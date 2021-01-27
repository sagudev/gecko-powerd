/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jit/ppc64/Assembler-ppc64.h"
#include "jit/FlushICache.h"

#include "mozilla/DebugOnly.h"

using mozilla::DebugOnly;

using namespace js;
using namespace js::jit;

ABIArgGenerator::ABIArgGenerator()
  : stackOffset_(0),
    usedGPRs_(0),
    usedFPRs_(0),
    current_()
{}

ABIArg
ABIArgGenerator::next(MIRType type)
{
    switch (type) {
      case MIRType::Int32:
      case MIRType::Int64:
      case MIRType::Pointer: {
        if (usedGPRs_ == 8)
            MOZ_CRASH("ABIArgGenerator overflow (GPR)");
        current_ = ABIArg(Register::FromCode((Register::Code)(usedGPRs_ + 3)));
        usedGPRs_++;
        break;
      }
      case MIRType::Float32:
      case MIRType::Double: {
        if (usedFPRs_ == 12)
            MOZ_CRASH("ABIArgGenerator overflow(FPRs)");

        current_ = ABIArg(FloatRegister::FromCode((Register::Code)(usedFPRs_ + 1)));
        usedGPRs_++;
        usedFPRs_++;
        break;
      }
      default:
        MOZ_CRASH("Unexpected argument type");
    }
    return current_;
}

uintptr_t
Assembler::GetPointer(uint8_t* instPtr)
{
    Instruction* inst = (Instruction*)instPtr;
    return Assembler::ExtractLoad64Value(inst);
}

static JitCode *
CodeFromJump(Instruction* jump)
{
    uint8_t* target = (uint8_t*)Assembler::ExtractLoad64Value(jump);
    return JitCode::FromExecutable(target);
}

void
Assembler::TraceJumpRelocations(JSTracer* trc, JitCode* code, CompactBufferReader& reader)
{
    while (reader.more()) {
        JitCode* child = CodeFromJump((Instruction*)(code->raw() + reader.readUnsigned()));
        TraceManuallyBarrieredEdge(trc, &child, "rel32");
    }
}

static void
TraceOneDataRelocation(JSTracer* trc, Instruction* inst)
{
    void* ptr = (void*)Assembler::ExtractLoad64Value(inst);
    void* prior = ptr;

    // All pointers on MIPS64 will have the top bits cleared. If those bits
    // are not cleared, this must be a Value.
    uintptr_t word = reinterpret_cast<uintptr_t>(ptr);
    if (word >> JSVAL_TAG_SHIFT) {
        Value v = Value::fromRawBits(word);
        TraceManuallyBarrieredEdge(trc, &v, "ion-masm-value");
        ptr = (void*)v.bitsAsPunboxPointer();
    } else {
        // No barrier needed since these are constants.
        TraceManuallyBarrieredGenericPointerEdge(trc, reinterpret_cast<gc::Cell**>(&ptr),
                                                     "ion-masm-ptr");
    }

    if (ptr != prior) {
        Assembler::UpdateLoad64Value(inst, uint64_t(ptr));
        FlushICache(inst, 6 * sizeof(uint32_t));
    }
}

/* static */ void
Assembler::TraceDataRelocations(JSTracer* trc, JitCode* code, CompactBufferReader& reader)
{
    while (reader.more()) {
        size_t offset = reader.readUnsigned();
        Instruction* inst = (Instruction*)(code->raw() + offset);
        TraceOneDataRelocation(trc, inst);
    }
}

void
Assembler::Bind(uint8_t* rawCode, const CodeLabel& label)
{
    if (label.patchAt().bound()) {

        auto mode = label.linkMode();
        intptr_t offset = label.patchAt().offset();
        intptr_t target = label.target().offset();

        if (mode == CodeLabel::RawPointer) {
            *reinterpret_cast<const void**>(rawCode + offset) = rawCode + target;
        } else {
            MOZ_ASSERT(mode == CodeLabel::MoveImmediate || mode == CodeLabel::JumpImmediate);
            Instruction* inst = (Instruction*) (rawCode + offset);
            Assembler::UpdateLoad64Value(inst, (uint64_t)(rawCode + target));
        }
    }
}

void
Assembler::bind(InstImm* inst, uintptr_t branch, uintptr_t target)
{
#if 0 // TODO: Assembler::bind()
    int64_t offset = target - branch;
    InstImm inst_bgezal = InstImm(op_regimm, r0, rt_bgezal, BOffImm16(0));
    InstImm inst_beq = InstImm(PPC_bc, r0, r0, BOffImm16(0));

    // If encoded offset is 4, then the jump must be short
    if (BOffImm16(inst[0]).decode() == 4) {
        MOZ_ASSERT(BOffImm16::IsInRange(offset));
        inst[0].setBOffImm16(BOffImm16(offset));
        inst[1].makeOp_nop();
        return;
    }

    // Generate the long jump for calls because return address has to be the
    // address after the reserved block.
    if (inst[0].encode() == inst_bgezal.encode()) {
        addLongJump(BufferOffset(branch));
        Assembler::WriteLoad64Instructions(inst, ScratchRegister, LabelBase::INVALID_OFFSET);
        inst[4] = InstReg(PPC_b | LinkB, ScratchRegister, r0).encode();
        // There is 1 nop after this.
        return;
    }

    if (BOffImm16::IsInRange(offset)) {
        // Don't skip trailing nops can improve performance
        // on Loongson3 platform.
        bool skipNops = (inst[0].encode() != inst_bgezal.encode() &&
                        inst[0].encode() != inst_beq.encode());

        inst[0].setBOffImm16(BOffImm16(offset));
        inst[1].makeOp_nop();

        if (skipNops) {
            inst[2] = InstImm(op_regimm, r0, rt_bgez, BOffImm16(5 * sizeof(uint32_t))).encode();
            // There are 4 nops after this
        }
        return;
    }

    if (inst[0].encode() == inst_beq.encode()) {
        // Handle long unconditional jump.
        addLongJump(BufferOffset(branch));
        Assembler::WriteLoad64Instructions(inst, ScratchRegister, LabelBase::INVALID_OFFSET);
        inst[4] = InstReg(op_special, ScratchRegister, r0, r0, ff_jr).encode();
        // There is 1 nop after this.
    } else {
        // Handle long conditional jump.
        inst[0] = invertBranch(inst[0], BOffImm16(7 * sizeof(uint32_t)));
        // No need for a "nop" here because we can clobber scratch.
        addLongJump(BufferOffset(branch + sizeof(uint32_t)));
        Assembler::WriteLoad64Instructions(&inst[1], ScratchRegister, LabelBase::INVALID_OFFSET);
        inst[5] = InstReg(op_special, ScratchRegister, r0, r0, ff_jr).encode();
        // There is 1 nop after this.
    }
#endif
}

#if 0
void
Assembler::bind(RepatchLabel* label)
{
    BufferOffset dest = nextOffset();
    if (label->used() && !oom()) {
        // If the label has a use, then change this use to refer to
        // the bound label;
        BufferOffset b(label->offset());
        InstImm* inst = (InstImm*)editSrc(b);
        InstImm inst_beq = InstImm(op_beq, r0, r0, BOffImm16(0));
        uint64_t offset = dest.getOffset() - label->offset();

        // If first instruction is lui, then this is a long jump.
        // If second instruction is lui, then this is a loop backedge.
        if (inst[0].extractOpcode() == (uint32_t(op_lui) >> OpcodeShift)) {
            // For unconditional long branches generated by ma_liPatchable,
            // such as under:
            //     jumpWithpatch
            addLongJump(BufferOffset(label->offset()));
        } else if (inst[1].extractOpcode() == (uint32_t(op_lui) >> OpcodeShift) ||
                   BOffImm16::IsInRange(offset))
        {
            // Handle code produced by:
            //     backedgeJump
            //     branchWithCode
            MOZ_ASSERT(BOffImm16::IsInRange(offset));
            MOZ_ASSERT(inst[0].extractOpcode() == (uint32_t(op_beq) >> OpcodeShift) ||
                       inst[0].extractOpcode() == (uint32_t(op_bne) >> OpcodeShift) ||
                       inst[0].extractOpcode() == (uint32_t(op_blez) >> OpcodeShift) ||
                       inst[0].extractOpcode() == (uint32_t(op_bgtz) >> OpcodeShift) ||
                       (inst[0].extractOpcode() == (uint32_t(op_regimm) >> OpcodeShift) &&
                       inst[0].extractRT() == (uint32_t(rt_bltz) >> RTShift)));
            inst[0].setBOffImm16(BOffImm16(offset));
        } else if (inst[0].encode() == inst_beq.encode()) {
            // Handle open long unconditional jumps created by
            // MacroAssemblerMIPSShared::ma_b(..., wasm::Trap, ...).
            // We need to add it to long jumps array here.
            // See MacroAssemblerMIPS64::branchWithCode().
            MOZ_ASSERT(inst[1].encode() == NopInst);
            MOZ_ASSERT(inst[2].encode() == NopInst);
            MOZ_ASSERT(inst[3].encode() == NopInst);
            MOZ_ASSERT(inst[4].encode() == NopInst);
            MOZ_ASSERT(inst[5].encode() == NopInst);
            addLongJump(BufferOffset(label->offset()));
            Assembler::WriteLoad64Instructions(inst, ScratchRegister, LabelBase::INVALID_OFFSET);
            inst[4] = InstReg(op_special, ScratchRegister, r0, r0, ff_jr).encode();
        } else {
            // Handle open long conditional jumps created by
            // MacroAssemblerMIPSShared::ma_b(..., wasm::Trap, ...).
            inst[0] = invertBranch(inst[0], BOffImm16(7 * sizeof(uint32_t)));
            // No need for a "nop" here because we can clobber scratch.
            // We need to add it to long jumps array here.
            // See MacroAssemblerMIPS64::branchWithCode().
            MOZ_ASSERT(inst[1].encode() == NopInst);
            MOZ_ASSERT(inst[2].encode() == NopInst);
            MOZ_ASSERT(inst[3].encode() == NopInst);
            MOZ_ASSERT(inst[4].encode() == NopInst);
            MOZ_ASSERT(inst[5].encode() == NopInst);
            MOZ_ASSERT(inst[6].encode() == NopInst);
            addLongJump(BufferOffset(label->offset() + sizeof(uint32_t)));
            Assembler::WriteLoad64Instructions(&inst[1], ScratchRegister, LabelBase::INVALID_OFFSET);
            inst[5] = InstReg(op_special, ScratchRegister, r0, r0, ff_jr).encode();
        }
    }
    label->bind(dest.getOffset());
}
#endif

void
Assembler::processCodeLabels(uint8_t* rawCode)
{
    for (const CodeLabel& label : codeLabels_) {
        Bind(rawCode, label);
    }
}

uint32_t
Assembler::PatchWrite_NearCallSize()
{
    // Load an address needs 5 instructions, mtctr, bctrl
    return (5 + 2) * sizeof(uint32_t);
}

void
Assembler::PatchWrite_NearCall(CodeLocationLabel start, CodeLocationLabel toCall)
{
    Instruction* inst = (Instruction*) start.raw();
    uint8_t* dest = toCall.raw();

    // Overwrite whatever instruction used to be here with a call.
    // Always use long jump for two reasons:
    // - Jump has to be the same size because of PatchWrite_NearCallSize.
    // - Return address has to be at the end of replaced block.
    // Short jump wouldn't be more efficient.
    Assembler::WriteLoad64Instructions(inst, ScratchRegister, (uint64_t)dest);
    // XXX: This is wrong!
    //inst[4] = InstReg(PPC_b | LinkB, ScratchRegister, r0, lr, ff_jalr);
    inst[5] = Instruction(PPC_nop);
    inst[6] = Instruction(PPC_nop);

    // Ensure everyone sees the code that was just written into memory.
    FlushICache(inst, PatchWrite_NearCallSize());
}

uint64_t
Assembler::ExtractLoad64Value(Instruction* inst0)
{
    InstImm* i0 = (InstImm*) inst0;
    InstImm* i1 = (InstImm*) i0->next();
    Instruction* i2 = (Instruction*) i1->next();
    InstImm* i3 = (InstImm*) i2->next();
    InstImm* i5 = (InstImm*) i3->next()->next();

    MOZ_ASSERT(i0->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i1->extractOpcode() == (uint32_t)PPC_ori);
    MOZ_ASSERT(i3->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i4->extractOpcode() == (uint32_t)PPC_ori);
    uint64_t value = (uint64_t(i0->extractImm16Value()) << 48) |
                     (uint64_t(i1->extractImm16Value()) << 32) |
                     (uint64_t(i3->extractImm16Value()) << 16) |
                     uint64_t(i5->extractImm16Value());
    return value;
}

void
Assembler::UpdateLoad64Value(Instruction* inst0, uint64_t value)
{
    InstImm* i0 = (InstImm*) inst0;
    InstImm* i1 = (InstImm*) i0->next();
    Instruction* i2 = (Instruction*) i1->next();
    InstImm* i3 = (InstImm*) i2->next();
    InstImm* i5 = (InstImm*) i3->next()->next();

    MOZ_ASSERT(i0->extractOpcode() == (uint32_t)op_oris);
    MOZ_ASSERT(i1->extractOpcode() == (uint32_t)PPC_ori);
    MOZ_ASSERT(i3->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i4->extractOpcode() == (uint32_t)PPC_ori);

    i0->setImm16(Imm16::Upper(Imm32(value >> 32)));
    i1->setImm16(Imm16::Lower(Imm32(value >> 32)));
    i3->setImm16(Imm16::Upper(Imm32(value)));
    i5->setImm16(Imm16::Lower(Imm32(value)));
}

void
Assembler::WriteLoad64Instructions(Instruction* inst0, Register reg, uint64_t value)
{
    Instruction* inst1 = inst0->next();
    Instruction* inst2 = inst1->next();
    Instruction* inst3 = inst2->next();
    Instruction* inst4 = inst3->next();

    *inst0 = InstImm(PPC_ori, reg, r0, BOffImm16(Imm16::Lower(Imm32(value >> 32)).encode()));
    *inst1 = InstImm(PPC_oris, reg, reg, BOffImm16(Imm16::Upper(Imm32(value >> 32)).encode()));
    // Manually construct 'sldi reg, reg, 32'
    *inst2 = InstImm(PPC_rldicr, reg, reg, BOffImm16((32 << 11)  | (31 << 4) | (1 << 2)));
    *inst3 = InstImm(PPC_oris, reg, reg, BOffImm16(Imm16::Upper(Imm32(value)).encode()));
    *inst4 = InstImm(PPC_ori, reg, reg, BOffImm16(Imm16::Lower(Imm32(value)).encode()));
}

void
Assembler::PatchDataWithValueCheck(CodeLocationLabel label, ImmPtr newValue,
                                   ImmPtr expectedValue)
{
    PatchDataWithValueCheck(label, PatchedImmPtr(newValue.value),
                            PatchedImmPtr(expectedValue.value));
}

void
Assembler::PatchDataWithValueCheck(CodeLocationLabel label, PatchedImmPtr newValue,
                                   PatchedImmPtr expectedValue)
{
    Instruction* inst = (Instruction*) label.raw();

    // Extract old Value
    DebugOnly<uint64_t> value = Assembler::ExtractLoad64Value(inst);
    MOZ_ASSERT(value == uint64_t(expectedValue.value));

    // Replace with new value
    Assembler::UpdateLoad64Value(inst, uint64_t(newValue.value));

    FlushICache(inst, 6 * sizeof(uint32_t));
}

void
Assembler::ToggleCall(CodeLocationLabel inst_, bool enabled)
{
    Instruction* inst = (Instruction*)inst_.raw();
    InstImm* i0 = (InstImm*) inst;
    InstImm* i1 = (InstImm*) i0->next();
    InstImm* i3 = (InstImm*) i1->next()->next();
    InstImm* i4 = (InstImm*) i3->next();
    Instruction* i5 = (Instruction*) i4->next();
    Instruction* i6 = (Instruction*) i5->next();

    MOZ_ASSERT(i0->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i1->extractOpcode() == (uint32_t)PPC_ori);
    MOZ_ASSERT(i3->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i4->extractOpcode() == (uint32_t)PPC_ori);

    if (enabled) {
        i5->makeOp_mtctr(ScratchRegister);
        i6->makeOp_bctr(LinkB);
    } else {
        *i5 = Instruction(PPC_nop);
        *i6 = Instruction(PPC_nop);
    }

    FlushICache(i4, sizeof(uint32_t));
}

bool Assembler::swapBuffer(wasm::Bytes& bytes) {
  // For now, specialize to the one use case. As long as wasm::Bytes is a
  // Vector, not a linked-list of chunks, there's not much we can do other
  // than copy.
  MOZ_ASSERT(bytes.empty());
  if (!bytes.resize(bytesNeeded())) {
    return false;
  }
  m_buffer.executableCopy(bytes.begin());
  return true;
}
