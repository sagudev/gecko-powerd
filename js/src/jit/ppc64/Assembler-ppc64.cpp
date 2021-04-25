/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jit/ppc64/Assembler-ppc64.h"
#include "jit/FlushICache.h"

#include "mozilla/DebugOnly.h"

#if DEBUG
#define spew(...) JitSpew(JitSpew_Codegen, __VA_ARGS__)
#else
#define spew(...)
#endif

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

bool Assembler::isGPR(Register r) const {
    return r.code() <  Registers::SPRStart;
}
bool Assembler::isSPR(Register r) const {
    return r.code() >= Registers::SPRStart;
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
    __asm__("trap\n");
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
    __asm__("trap\n");
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

void
Assembler::bind(Label* label)
{
    __asm__("trap\n");
#if 0
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
            // MacroAssemblerMIPSShared::ma_bc(..., wasm::Trap, ...).
            // We need to add it to long jumps array here.
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
            // MacroAssemblerMIPSShared::ma_bc(..., wasm::Trap, ...).
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
#endif
}

void
Assembler::bind(Label* label, BufferOffset boff)
{
    __asm__("trap\n");
    spew(".set Llabel %p", label);
    // If our caller didn't give us an explicit target to bind to
    // then we want to bind to the location of the next instruction
    BufferOffset dest = boff.assigned() ? boff : nextOffset();
    if (label->used()) {
        int32_t next;

        // A used label holds a link to branch that uses it.
        BufferOffset b(label);
        do {
            // Even a 0 offset may be invalid if we're out of memory.
            if (oom()) {
                return;
            }

            Instruction* inst = editSrc(b);

            // Second word holds a pointer to the next branch in label's chain.
            next = inst[1].encode();
            bind(reinterpret_cast<InstImm*>(inst), b.getOffset(), dest.getOffset());

            b = BufferOffset(next);
        } while (next != LabelBase::INVALID_OFFSET);
    }
    label->bind(dest.getOffset());
}
void
Assembler::retarget(Label* label, Label* target)
{
    spew("retarget %p -> %p", label, target);
    if (label->used() && !oom()) {
        if (target->bound()) {
            bind(label, BufferOffset(target));
        } else if (target->used()) {
            // The target is not bound but used. Prepend label's branch list
            // onto target's.
            int32_t next;
            BufferOffset labelBranchOffset(label);

            // Find the head of the use chain for label.
            do {
                Instruction* inst = editSrc(labelBranchOffset);

                // Second word holds a pointer to the next branch in chain.
                next = inst[1].encode();
                labelBranchOffset = BufferOffset(next);
            } while (next != LabelBase::INVALID_OFFSET);

            // Then patch the head of label's use chain to the tail of
            // target's use chain, prepending the entire use chain of target.
            Instruction* inst = editSrc(labelBranchOffset);
            int32_t prev = target->offset();
            target->use(label->offset());
            inst[1].setData(prev);
        } else {
            // The target is unbound and unused.  We can just take the head of
            // the list hanging off of label, and dump that into target.
            target->use(label->offset());
        }
    }
    label->reset();
}

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
    __asm__("trap\n");
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

void
Assembler::PatchWrite_Imm32(CodeLocationLabel label, Imm32 imm)
{
    uint32_t *l = (uint32_t *)label.raw();

    *(l - 1) = imm.value;
}

void
Assembler::UpdateLisOriValue(Instruction *inst0, Instruction *inst1,
                             uint32_t value)
{
    MOZ_ASSERT(inst0->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(inst1->extractOpcode() == (uint32_t)PPC_ori);

    reinterpret_cast<InstImm*>(inst0)->setImm16(value >> 16);
    reinterpret_cast<InstImm*>(inst1)->setImm16(value & 0xffff);
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
    MOZ_ASSERT(i5->extractOpcode() == (uint32_t)PPC_ori);
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

    MOZ_ASSERT(i0->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i1->extractOpcode() == (uint32_t)PPC_ori);
    MOZ_ASSERT(i3->extractOpcode() == (uint32_t)PPC_oris);
    MOZ_ASSERT(i5->extractOpcode() == (uint32_t)PPC_ori);

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

    *inst0 = InstImm(PPC_ori, reg, r0, Imm16::Lower(Imm32(value >> 32)).encode());
    *inst1 = InstImm(PPC_oris, reg, reg, Imm16::Upper(Imm32(value >> 32)).encode());
    // Manually construct 'sldi reg, reg, 32'
    *inst2 = InstImm(PPC_rldicr, reg, reg, Imm16((32 << 11)  | (31 << 4) | (1 << 2)));
    *inst3 = InstImm(PPC_oris, reg, reg, Imm16::Upper(Imm32(value)).encode());
    *inst4 = InstImm(PPC_ori, reg, reg, Imm16::Lower(Imm32(value)).encode());
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

void
Assembler::ToggleToJmp(CodeLocationLabel inst_)
{
    InstImm* inst = (InstImm*)inst_.raw();

    __asm__("trap\n");
    MOZ_ASSERT(inst->extractOpcode() == PPC_andis);
#if 0 // TODO
    // We converted beq to andi, so now we restore it.
    inst->setOpcode(op_beq);
    *inst = InstImm(PPC_bc, 
#endif
}

void
Assembler::ToggleToCmp(CodeLocationLabel inst_)
{
    InstImm* inst = (InstImm*)inst_.raw();

    __asm__("trap\n");
    // toggledJump is allways used for short jumps.
    MOZ_ASSERT(inst->extractOpcode() == PPC_bc);
#if 0 // TODO
    // Replace "beq $zero, $zero, offset" with "andi $zero, $zero, offset"
    inst->setOpcode(op_andi);
#endif
}

Assembler::Condition
Assembler::InvertCondition( Condition cond)
{
    switch (cond) {
        case Equal:
            return NotEqual;
        case NotEqual:
            return Equal;
        case LessThan:
            return GreaterThanOrEqual;
        case LessThanOrEqual:
            return GreaterThan;
        case GreaterThan:
            return LessThanOrEqual;
        case GreaterThanOrEqual:
            return LessThan;
        case Above:
            return BelowOrEqual;
        case AboveOrEqual:
            return Below;
        case Below:
            return AboveOrEqual;
        case BelowOrEqual:
            return Above;
        default:
            MOZ_CRASH("unexpected condition");
    }
}

Assembler::DoubleCondition
Assembler::InvertCondition( DoubleCondition cond)
{
    switch (cond) {
        case DoubleOrdered:
            return DoubleUnordered;
        case DoubleEqual:
            return DoubleNotEqualOrUnordered;
        case DoubleNotEqual:
            return DoubleEqualOrUnordered;
        case DoubleGreaterThan:
            return DoubleLessThanOrEqualOrUnordered;
        case DoubleGreaterThanOrEqual:
            return DoubleLessThanOrUnordered;
        case DoubleLessThan:
            return DoubleGreaterThanOrEqualOrUnordered;
        case DoubleLessThanOrEqual:
            return DoubleGreaterThanOrUnordered;
        case DoubleUnordered:
            return DoubleOrdered;
        case DoubleEqualOrUnordered:
            return DoubleNotEqual;
        case DoubleNotEqualOrUnordered:
            return DoubleEqual;
        case DoubleGreaterThanOrUnordered:
            return DoubleLessThanOrEqual;
        case DoubleGreaterThanOrEqualOrUnordered:
            return DoubleLessThan;
        case DoubleLessThanOrUnordered:
            return DoubleGreaterThanOrEqual;
        case DoubleLessThanOrEqualOrUnordered:
            return DoubleGreaterThan;
        default:
            MOZ_CRASH("unexpected condition");
    }
}

bool
Assembler::swapBuffer(wasm::Bytes& bytes) {
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

void
Assembler::copyJumpRelocationTable(uint8_t* dest)
{
    if (jumpRelocations_.length()) {
        memcpy(dest, jumpRelocations_.buffer(), jumpRelocations_.length());
    }
}

void
Assembler::copyDataRelocationTable(uint8_t* dest)
{
    if (dataRelocations_.length()) {
        memcpy(dest, dataRelocations_.buffer(), dataRelocations_.length());
    }
}

void
Assembler::finish()
{
    MOZ_ASSERT(!isFinished);
    isFinished = true;
}

void
Assembler::executableCopy(uint8_t* buffer)
{
    MOZ_ASSERT(isFinished);
    m_buffer.executableCopy(buffer);
}

bool
Assembler::appendRawCode(unsigned char const *bytes, unsigned long length)
{
    return m_buffer.appendRawCode(bytes, length);
}

bool Assembler::oom() const {
  return AssemblerShared::oom() || m_buffer.oom() || jumpRelocations_.oom() ||
         dataRelocations_.oom();
}

BufferOffset Assembler::haltingAlign(int align)
{
    BufferOffset ret;
    MOZ_ASSERT(m_buffer.isAligned(4));
    if (align == 8) {
        if (!m_buffer.isAligned(align)) {
            BufferOffset tmp = xs_trap();
            if (!ret.assigned()) {
                ret = tmp;
            }
        }
    } else {
        MOZ_ASSERT((align& (align- 1)) == 0);
        while (size() & (align- 1)) {
            BufferOffset tmp = xs_trap();
            if (!ret.assigned()) {
                ret = tmp;
            }
        }
    }
    return ret;
}

BufferOffset Assembler::nopAlign(int align)
{
    BufferOffset ret;
    MOZ_ASSERT(m_buffer.isAligned(4));
    if (align == 8) {
        if (!m_buffer.isAligned(align)) {
            BufferOffset tmp = as_nop();
            if (!ret.assigned()) {
                ret = tmp;
            }
        }
    } else {
        MOZ_ASSERT((align& (align- 1)) == 0);
        while (size() & (align- 1)) {
            BufferOffset tmp = as_nop();
            if (!ret.assigned()) {
                ret = tmp;
            }
        }
    }
    return ret;
}

size_t
Assembler::size() const
{
    return m_buffer.size();
}

size_t
Assembler::dataRelocationTableBytes() const
{
    return dataRelocations_.length();
}

size_t
Assembler::jumpRelocationTableBytes() const
{
    return jumpRelocations_.length();
}

size_t
Assembler::bytesNeeded() const
{
    return m_buffer.size() +
        jumpRelocations_.length() +
        dataRelocations_.length() +
        preBarriers_.length();
}

BufferOffset
Assembler::writeInst(uint32_t x, uint32_t *dest)
{
    if (dest == nullptr)
        return m_buffer.putInt(x);

    *dest = x;
    return BufferOffset();
}

BufferOffset Assembler::as_nop()
{
    spew("nop");
    return writeInst(PPC_nop);
}

BufferOffset Assembler::as_lwsync()
{
    spew("lwsync");
    return writeInst(PPC_lwsync);
}

BufferOffset Assembler::as_sync()
{
    spew("sync");
    return writeInst(PPC_sync);
}

// Branch and jump instructions.
BufferOffset Assembler::as_b(JOffImm26 off, BranchAddressType bat, LinkBit lb)
{
    return as_b(off.encode(), bat, lb);
}

BufferOffset Assembler::as_b(int32_t off, BranchAddressType bat, LinkBit lb)
{
    spew("b%s%s\t%x\n", bat == AbsoluteBranch ? "a" : "", lb ? "l" : "", off);
    return writeInst(PPC_b | off | bat | lb);
}

BufferOffset Assembler::as_blr(LinkBit lb)
{
    spew("blr%s", lb ? "l" : "");
    return writeInst(PPC_blr | lb);
}

BufferOffset Assembler::as_bctr(LinkBit lb)
{
    spew("bctr%s", lb ? "l" : "");
    return writeInst(PPC_blr | lb);
}
    
static uint32_t makeOpMask(Assembler::DoubleCondition cond, CRegisterID cr)
{
}

static uint32_t makeOpMask(Assembler::Condition cond, CRegisterID cr)
{
}

// Conditional branches.
BufferOffset Assembler::as_bc(BOffImm16 off, Condition cond, CRegisterID cr,
        LikelyBit lkb, LinkBit lb)
{
    return as_bc(off.encode(), cond, cr, lkb, lb);
}

BufferOffset Assembler::as_bc(int16_t off, Condition cond, CRegisterID cr,
        LikelyBit lkb, LinkBit lb)
{
    return as_bc(off, makeOpMask(cond, cr), lkb, lb);
}

BufferOffset Assembler::as_bc(BOffImm16 off, DoubleCondition cond,
        CRegisterID cr, LikelyBit lkb, LinkBit lb)
{
    return as_bc(off.encode(), cond, cr, lkb, lb);
}

BufferOffset Assembler::as_bc(int16_t off, DoubleCondition cond, CRegisterID cr,
        LikelyBit lkb, LinkBit lb)
{
    return as_bc(off, makeOpMask(cond, cr), lkb, lb);
}

BufferOffset Assembler::as_bcctr(Condition cond, CRegisterID cr, LikelyBit lkb,
        LinkBit lb)
{
    return as_bcctr(makeOpMask(cond, cr), lkb, lb);
}

BufferOffset Assembler::as_bcctr(DoubleCondition cond, CRegisterID cr,
        LikelyBit lkb, LinkBit lb)
{
    return as_bcctr(makeOpMask(cond, cr), lkb, lb);
}
    
static uint32_t makeOpMask(uint32_t op)
{
    return ((op & 0x0f) << 21) | ((op & 0xf0) << 12);
}

BufferOffset Assembler::as_bc(int16_t off, uint32_t op, LikelyBit lkb, LinkBit lb)
{
    spew("bc");
    return writeInst(Instruction(PPC_bc | makeOpMask(op) | lkb << 21 | off | lb).encode());
}

BufferOffset Assembler::as_bcctr(uint32_t op, LikelyBit lkb, LinkBit lb)
{
    spew("bcctr");
    return writeInst(PPC_bcctr | makeOpMask(op) | lkb << 21 | lb);
}

// SPR operations.
BufferOffset Assembler::as_mtspr(SPRegisterID spr, Register ra)
{
    spew("mtspr\t%d,%3s", spr, ra.name());
    return writeInst(PPC_mfspr | ra.code() << 21 | PPC_SPR(spr) << 11);
}
BufferOffset Assembler::as_mfspr(Register rd, SPRegisterID spr)
{
    spew("mfspr\t%3s,%d", rd.name(), spr);
    return writeInst(PPC_mfspr | rd.code() << 21 | PPC_SPR(spr) << 11);
}

// CR operations.
#define DEF_CRCR(op) \
        BufferOffset Assembler::as_##op(uint8_t t, uint8_t a, uint8_t b) { \
            spew(#op"\t%d,%d,%d", t, a, b); \
            return writeInst(PPC_##op | t << 21 | a << 16 | b << 11); }
DEF_CRCR(crand)
DEF_CRCR(crandc)
DEF_CRCR(cror)
DEF_CRCR(crorc)
DEF_CRCR(crxor)
#undef DEF_CRCR

BufferOffset Assembler::as_mtcrf(uint32_t mask, Register rs)
{
    spew("mtcrf %d,%3s", mask, rs.name());
    return writeInst(PPC_mtcrf | rs.code() << 21 | mask << 11);
}

BufferOffset Assembler::as_mfcr(Register rd)
{
    spew("mfcr %3s", rd.name());
    return writeInst(PPC_mfcr | rd.code() << 21);
}

BufferOffset Assembler::as_mfocrf(Register rd, CRegisterID crfs)
{
    spew("mfocrf %3s,cr%d", rd.name(), crfs);
    return writeInst(PPC_mfocrf | rd.code() << 21 | crfs << 12);
}

// XXX: This should now be mcrxrx
BufferOffset Assembler::as_mcrxr(CRegisterID crt, Register temp)
{
    spew("mcrxr\tcr%d", crt);
    return writeInst(PPC_mcrxr | crt << 23);
}

// GPR operations and load-stores.
BufferOffset Assembler::as_neg(Register rd, Register rs)
{
    spew("neg %3s,%3s", rd.name(), rs.name());
    return writeInst(InstReg(PPC_neg, rd, rs, r0).encode());
}

BufferOffset Assembler::as_cmpd(CRegisterID cr, Register ra, Register rb)
{
    spew("cmpd\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmpd | cr << 23 | ra.code() << 16 | rb.code() << 11);
}

BufferOffset Assembler::as_cmpdi(CRegisterID cr, Register ra, int16_t im)
{
    spew("cmpdi\tcr%d,%3s,%d", cr, ra.name(), im);
    return writeInst(PPC_cmpdi | cr << 23 | ra.code() << 16 | im);
}

BufferOffset Assembler::as_cmpld(CRegisterID cr, Register ra, Register rb)
{
    spew("cmpld\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmpld | cr << 23 | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmpldi(CRegisterID cr, Register ra, int16_t im)
{
    spew("cmpldi\tcr%d,%3s,%d", cr, ra.name(), im);
    return writeInst(PPC_cmpldi | cr << 23 | ra.code() << 16 | im);
}
BufferOffset Assembler::as_cmpw(CRegisterID cr, Register ra, Register rb)
{
    spew("cmpw\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmpw | cr << 23 | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmpwi(CRegisterID cr, Register ra, int16_t im)
{
    spew("cmpwi\tcr%d,%3s,%d", cr, ra.name(), im);
    return writeInst(PPC_cmpwi | cr << 23 | ra.code() << 16 | im);
}
BufferOffset Assembler::as_cmplw(CRegisterID cr, Register ra, Register rb)
{
    spew("cmplw\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmplw | cr << 23 | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmplwi(CRegisterID cr, Register ra, int16_t im)
{
    spew("cmplwi\tcr%d,%3s,%d", cr, ra.name(), im);
    return writeInst(PPC_cmplwi | cr << 23 | ra.code() << 16 | im);
}
BufferOffset Assembler::as_cmpd(Register ra, Register rb)
{
    spew("cmpd\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmpd | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmpdi(Register ra, int16_t im)
{
    spew("cmpdi\t%3s,%d", ra.name(), im);
    return writeInst(PPC_cmpdi | ra.code() << 16 | im);
}
BufferOffset Assembler::as_cmpld(Register ra, Register rb)
{
    spew("cmpld\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmpld | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmpldi(Register ra, int16_t im)
{
    spew("cmpdi\t%3s,%d", ra.name(), im);
    return writeInst(PPC_cmpdi | ra.code() << 16 | im);
}
BufferOffset Assembler::as_cmpw(Register ra, Register rb)
{
    spew("cmpw\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmpw | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmpwi(Register ra, int16_t im)
{
    spew("cmpwi\t%3s,%d", ra.name(), im);
    return writeInst(PPC_cmpwi | ra.code() << 16 | im);
}

BufferOffset Assembler::as_cmplw(Register ra, Register rb)
{
    spew("cmplw\t%3s,%3s", ra.name(), rb.name());
    return writeInst(PPC_cmplw | ra.code() << 16 | rb.code() << 11);
}
BufferOffset Assembler::as_cmplwi(Register ra, int16_t im)
{
    spew("cmplwi\t%3s,%d", ra.name(), im);
    return writeInst(PPC_cmplwi | ra.code() << 16 | im);
}

static uint32_t
AForm(uint32_t op, FloatRegister frt, FloatRegister fra, FloatRegister frb,
        FloatRegister frc, bool rc)
{
    return (op | frt.code() << 21 | fra.code() << 16 | frb.code() << 11 |
            frc.code() << 6 | rc);
}

static uint32_t
XForm(uint32_t op, FloatRegister frt, FloatRegister fra, FloatRegister frb, bool rc)
{
    return (op | frt.code() << 21 | fra.code() << 16 | frb.code() << 11 | rc);
}

static uint32_t
XForm(uint32_t op, FloatRegister frt, Register ra, Register rb, bool rc)
{
    return (op | frt.code() << 21 | ra.code() << 16 | rb.code() << 11 | rc);
}

static uint32_t
DForm(uint32_t op, FloatRegister frt, Register ra, int16_t imm)
{
    return (op | frt.code() << 21 | ra.code() << 16 | imm);
}

#define DEF_XFORM(op)   \
    BufferOffset Assembler::as_##op(Register rd, Register ra, Register rb) { \
        spew(#op "\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(InstReg(PPC_##op, rd, ra, rb).encode()); }

#define DEF_XFORM_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register rd, Register ra, Register rb) {\
        spew(#op ".\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(InstReg(PPC_##op, rd, ra, rb).encode() | 0x1); }

#define DEF_XFORMS(op)   \
    BufferOffset Assembler::as_##op(Register rd, Register ra, Register rb) { \
        spew(#op "\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(InstReg(PPC_##op, ra, rd, rb).encode()); }

#define DEF_XFORMS_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register rd, Register ra, Register rb) {\
        spew(#op ".\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(InstReg(PPC_##op, ra, rd, rb).encode() | 0x1); }

#define DEF_AFORM_C(op)   \
    BufferOffset Assembler::as_##op(FloatRegister rd, FloatRegister ra, FloatRegister rc) { \
        spew(#op "\t%3s,%3s,%3s", rd.name(), ra.name(), rc.name()); \
        return writeInst(AForm(PPC_##op, rd, ra, f0, rc, false)); }

#define DEF_AFORM_C_RC(op) \
    BufferOffset Assembler::as_##op##_rc(FloatRegister rd, FloatRegister ra, FloatRegister rc) {\
        spew(#op ".\t%3s,%3s,%3s", rd.name(), ra.name(), rc.name()); \
        return writeInst(AForm(PPC_##op, rd, ra, f0, rc, true)); }

#define DEF_AFORM_B(op)   \
    BufferOffset Assembler::as_##op(FloatRegister rd, FloatRegister ra, FloatRegister rb) { \
        spew(#op "\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(AForm(PPC_##op, rd, ra, rb, f0, false)); }

#define DEF_AFORM_B_RC(op) \
    BufferOffset Assembler::as_##op##_rc(FloatRegister rd, FloatRegister ra, FloatRegister rb) {\
        spew(#op ".\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(AForm(PPC_##op, rd, ra, rb, f0, true)); }

#define DEF_AFORM(op)   \
    BufferOffset Assembler::as_##op(FloatRegister rd, FloatRegister ra, FloatRegister rc, FloatRegister rb) { \
        spew(#op "\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(AForm(PPC_##op, rd, ra, rb, rc, false)); }

#define DEF_AFORM_RC(op) \
    BufferOffset Assembler::as_##op##_rc(FloatRegister rd, FloatRegister ra, FloatRegister rc, FloatRegister rb) {\
        spew(#op ".\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(AForm(PPC_##op, rd, ra, rb, rc, true)); }

#define DEF_XFORMS_I(op)   \
    BufferOffset Assembler::as_##op(Register rd, Register ra, uint8_t sh) { \
        spew(#op "\t%3s,%3s,%d", rd.name(), ra.name(), sh); \
        return writeInst(PPC_##op | ra.code() << 21 | rd.code() << 16 | sh << 11); }

#define DEF_XFORMS_I_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register rd, Register ra, uint8_t sh) {\
        spew(#op ".\t%3s,%3s,%d", rd.name(), ra.name(), sh); \
        return writeInst(PPC_##op | ra.code() << 21 | rd.code() << 16 | sh << 11 | 0x1); }

#define DEF_XFORM2(op)   \
    BufferOffset Assembler::as_##op(Register rd, Register ra) { \
        spew(#op "\t%3s,%3s", rd.name(), ra.name()); \
        return writeInst(InstReg(PPC_##op, rd, ra, r0).encode()); }

#define DEF_XFORM2_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register rd, Register ra) {\
        spew(#op ".\t%3s,%3s", rd.name(), ra.name()); \
        return writeInst(InstReg(PPC_##op, rd, ra, r0).encode() | 0x1); }

#define DEF_XFORM2_F(op) \
    BufferOffset Assembler::as_##op(FloatRegister rd, FloatRegister ra) { \
        spew(#op "\t%3s,%3s", rd.name(), ra.name()); \
        return writeInst(XForm(PPC_##op, rd, f0, ra, false)); }

#define DEF_XFORM2_F_RC(op) \
    BufferOffset Assembler::as_##op##_rc(FloatRegister rd, FloatRegister ra) {\
        spew(#op ".\t%3s,%3s", rd.name(), ra.name()); \
        return writeInst(XForm(PPC_##op, rd, f0, ra, true)); }

#define DEF_XFORM2S(op)   \
    BufferOffset Assembler::as_##op(Register rd, Register ra) { \
        spew(#op "\t%3s,%3s", rd.name(), ra.name()); \
        return writeInst(InstReg(PPC_##op, ra, rd, r0).encode()); }

#define DEF_XFORM2S_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register rd, Register ra) {\
        spew(#op ".\t%3s,%3s", rd.name(), ra.name()); \
        return writeInst(InstReg(PPC_##op, ra, rd, r0).encode() | 0x1); }

#define DEF_DFORM(op)   \
    BufferOffset Assembler::as_##op(Register ra, Register rs, int16_t im) { \
        spew(#op "\t%3s,%3s,%d", ra.name(), rs.name(), im); \
        return writeInst(InstImm(PPC_##op, rs, ra, im).encode()); }

#define DEF_DFORM_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register ra, Register rs, uint16_t im) {\
        spew(#op ".\t%3s,%3s,%d", ra.name(), rs.name(), im); \
        return writeInst(InstImm(PPC_##op, rs, ra, im).encode() | 0x1); }

#define DEF_DFORMS(op)   \
    BufferOffset Assembler::as_##op(Register ra, Register rs, uint16_t im) { \
        spew(#op "\t%3s,%3s,%d", ra.name(), rs.name(), im); \
        return writeInst(InstImm(PPC_##op, ra, rs, im).encode()); }

#define DEF_DFORMS_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register ra, Register rs, uint16_t im) {\
        spew(#op ".\t%3s,%3s,%d", ra.name(), rs.name(), im); \
        return writeInst(InstImm(PPC_##op, ra, rs, im).encode() | 0x1); }

#define DEF_DFORM_F(op)   \
    BufferOffset Assembler::as_##op(FloatRegister rt, Register ra, int16_t im) { \
        spew(#op "\t%3s,%3s,%d", rt.name(), ra.name(), im); \
        return writeInst(DForm(PPC_##op, rt, ra, im)); }

#define DEF_MFORM(op) \
    BufferOffset Assembler::as_##op(Register ra, Register rs, Register rb, uint8_t mb, uint8_t me) {\
        spew(#op "\t%3s,%3s,%3s,%d,%d", ra.name(), rs.name(), rb.name(), mb, me); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | rb.code() << 11 | mb << 6 | me << 1); }

#define DEF_MFORM_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register ra, Register rs, Register rb, uint8_t mb, uint8_t me) {\
        spew(#op ".\t%3s,%3s,%3s,%d,%d", ra.name(), rs.name(), rb.name(), mb, me); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | rb.code() << 11 | mb << 6 | me << 1 | 1); }

#define DEF_MFORM_I(op) \
    BufferOffset Assembler::as_##op(Register ra, Register rs, uint8_t sh, uint8_t mb, uint8_t me) {\
        spew(#op "\t%3s,%3s,%d,%d,%d", ra.name(), rs.name(), sh, mb, me); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | sh << 11 | mb << 6 | me << 1); }

#define DEF_MFORM_I_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register ra, Register rs, uint8_t sh, uint8_t mb, uint8_t me) {\
        spew(#op ".\t%3s,%3s,%d,%d,%d", ra.name(), rs.name(), sh, mb, me); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | sh << 11 | mb << 6 | me << 1 | 1); }

#define DEF_MDSFORM(op) \
    BufferOffset Assembler::as_##op(Register ra, Register rs, Register rb, uint8_t mb) {\
        spew(#op "\t%3s,%3s,%3s,%d", ra.name(), rs.name(), rb.name(), mb); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | rb.code() << 11 | mb << 6); }

#define DEF_MDSFORM_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register ra, Register rs, Register rb, uint8_t mb) {\
        spew(#op ".\t%3s,%3s,%3s,%d", ra.name(), rs.name(), rb.name(), mb); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | rb.code() << 11 | mb << 6); }

#define DEF_MDFORM(op) \
    BufferOffset Assembler::as_##op(Register ra, Register rs, uint8_t sh, uint8_t mb) {\
        spew(#op "\t%3s,%3s,%d,%d", ra.name(), rs.name(), sh, mb); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | sh << 11 | mb << 6); }

#define DEF_MDFORM_RC(op) \
    BufferOffset Assembler::as_##op##_rc(Register ra, Register rs, uint8_t sh, uint8_t mb) {\
        spew(#op ".\t%3s,%3s,%d,%d", ra.name(), rs.name(), sh, mb); \
        return writeInst(PPC_##op | rs.code() << 21 | ra.code() << 16 | sh << 11 | mb << 6); }



DEF_MFORM(rlwnm)
DEF_MFORM_I(rlwinm)
DEF_MFORM_I_RC(rlwinm)
DEF_MFORM_I(rlwimi)
DEF_XFORMS_I(srawi)

DEF_MDSFORM(rldcl)
//DEF_MDSFORM_RC(rldcl)
//DEF_MDSFORM(rldcr)
DEF_MDFORM(rldicl)
DEF_MDFORM_RC(rldicl)
DEF_MDFORM(rldicr)
DEF_MDFORM_RC(rldicr)
DEF_MDFORM(rldimi)
//DEF_MDFORM_RC(rldimi)
BufferOffset Assembler::as_sradi(Register rd, Register rs, int sh)
{
    spew("sradi\t%3s,%3s,%d", rd.name(), rs.name(), sh);
    return writeInst(PPC_sradi | rd.code() << 16 | rs.code() << 21 |
            (sh & 0x1f) << 11 | (sh & 0x20) >> 4);
}

#define DEF_ALU2(op)    DEF_XFORM(op) DEF_XFORM_RC(op)

DEF_ALU2(add)
DEF_ALU2(addc)
DEF_ALU2(adde)
DEF_ALU2(addo)
DEF_ALU2(subf)
DEF_ALU2(subfc)
DEF_ALU2(subfe)
DEF_ALU2(subfo)
DEF_ALU2(divd)
DEF_ALU2(divdo)
DEF_ALU2(divdu)
DEF_ALU2(divduo)
DEF_ALU2(divw)
DEF_ALU2(divwo)
DEF_ALU2(divwu)
DEF_ALU2(divwuo)
DEF_ALU2(mulld)
DEF_ALU2(mulhd)
DEF_ALU2(mulhdu)
DEF_ALU2(mulldo)
DEF_ALU2(mullw)
DEF_ALU2(mulhw)
DEF_ALU2(mulhwu)
DEF_ALU2(mullwo)
DEF_ALU2(eqv) // NB: Implemented differently.
#undef DEF_ALU2

#define DEF_ALUI(op)    \
    BufferOffset Assembler::as_##op(Register rd, Register ra, int16_t im) { \
        spew(#op "\t%3s,%3s,%d", rd.name(), ra.name(), im); \
        return writeInst(InstImm(PPC_##op, rd, ra, im).encode()); }\
    BufferOffset Assembler::as_##op##_rc(Register rd, Register ra, int16_t im) { \
        spew(#op ".\t%3s,%3s,%d", rd.name(), ra.name(), im); \
        return writeInst(InstImm(PPC_##op, rd, ra, im).encode() | 0x1); }
// mscdfr0
BufferOffset Assembler::as_addi(Register rd, Register ra, int16_t im, bool actually_li) {
#if DEBUG
    if (actually_li) {
        spew("li\t%3s,%d", rd.name(), im);
    } else {
        MOZ_ASSERT(ra != r0); // Because that would be li
        spew("addi\t%3s,%3s,%d", rd.name(), ra.name(), im);
    }
#endif
    return writeInst(InstImm(PPC_addi, rd, ra, im).encode());
}
BufferOffset Assembler::as_addi_rc(Register rd, Register ra, int16_t im, bool actually_li) {
#if DEBUG
    if (actually_li) {
        spew("li.\t%3s,%d", rd.name(), im);
    } else {
        MOZ_ASSERT(ra != r0); // Because that would be li
        spew("addi.\t%3s,%3s,%d", rd.name(), ra.name(), im);
    }
#endif
    return writeInst(InstImm(PPC_addi, rd, ra, im).encode() | 0x01);
}
BufferOffset Assembler::as_addis(Register rd, Register ra, int16_t im, bool actually_lis) {
#if DEBUG
    if (actually_lis) {
        spew("lis\t%3s,%d", rd.name(), im);
    } else {
        MOZ_ASSERT(ra != r0); // Because that would be lis
        spew("addis\t%3s,%3s,%d", rd.name(), ra.name(), im);
    }
#endif
    return writeInst(InstImm(PPC_addis, rd, ra, im).encode());
}
BufferOffset Assembler::as_addis_rc(Register rd, Register ra, int16_t im, bool actually_lis) {
#if DEBUG
    if (actually_lis) {
        spew("lis.\t%3s,%d", rd.name(), im);
    } else {
        MOZ_ASSERT(ra != r0); // Because that would be lis
        spew("addis.\t%3s,%3s,%d", rd.name(), ra.name(), im);
    }
#endif
    return writeInst(InstImm(PPC_addis, rd, ra, im).encode() | 0x01);
}

DEF_ALUI(addic)
// NB: mulli is usually strength-reduced, since it can take up to five
// cycles in the worst case.
DEF_ALUI(mulli)
DEF_ALUI(subfic)
#undef DEF_ALUI

#define DEF_ALUE(op) DEF_XFORM2(op) DEF_XFORM2_RC(op)
DEF_ALUE(addme)
DEF_ALUE(addze)
DEF_ALUE(subfze)
#undef DEF_ALUE

#define DEF_ALUE(op) DEF_XFORM2S(op) DEF_XFORM2S_RC(op)
DEF_ALUE(cntlzw) // NB: In this case, rd = ra and ra = rs, but no biggie here.
DEF_ALUE(cntlzd) // NB: In this case, rd = ra and ra = rs, but no biggie here.
DEF_ALUE(cnttzd) // NB: In this case, rd = ra and ra = rs, but no biggie here.
DEF_ALUE(cnttzw) // NB: In this case, rd = ra and ra = rs, but no biggie here.
#undef DEF_ALUE

DEF_XFORM2S(popcntd)
DEF_XFORM2S(popcntw)

#define DEF_BITALU2(op) DEF_XFORMS(op) DEF_XFORMS_RC(op)
DEF_BITALU2(andc)
DEF_BITALU2(nand)
DEF_BITALU2(nor)
DEF_BITALU2(slw)
DEF_BITALU2(srw)
DEF_BITALU2(sraw)
DEF_BITALU2(sld)
DEF_BITALU2(srd)
DEF_BITALU2(srad)
DEF_BITALU2(and)
DEF_BITALU2(or)
DEF_BITALU2(xor)
#undef DEF_BITALU2

#define DEF_BITALUI(op) DEF_DFORMS(op) DEF_DFORMS_RC(op)
DEF_BITALUI(ori)
DEF_BITALUI(oris)
DEF_BITALUI(xori)
DEF_BITALUI(xoris)
#undef DEF_BITALUI
DEF_DFORMS_RC(andi)
DEF_DFORMS_RC(andis)
        
#define DEF_ALUEXT(op) DEF_XFORM2S(op) DEF_XFORM2S_RC(op)
DEF_ALUEXT(extsb)
DEF_ALUEXT(extsh)
DEF_ALUEXT(extsw)
#undef DEF_ALUEXT

#define DEF_MEMd(op) DEF_DFORM(op)
DEF_MEMd(lbz)
DEF_MEMd(lha)
DEF_MEMd(lhz)
DEF_MEMd(lwz)
DEF_MEMd(ld)

DEF_MEMd(stb)
DEF_MEMd(stw)
DEF_MEMd(stwu)
DEF_MEMd(sth)
DEF_MEMd(std)
DEF_MEMd(stdu)
#undef DEF_MEMd

#define DEF_MEMx(op) DEF_XFORM(op)
DEF_MEMx(lbzx)
DEF_MEMx(lhax)
DEF_MEMx(lhzx)
DEF_MEMx(lhbrx)
DEF_MEMx(lwzx)
DEF_MEMx(lwbrx)
DEF_MEMx(lwarx)
DEF_MEMx(ldx)
DEF_MEMx(ldarx)

DEF_MEMx(stbx)
DEF_MEMx(stwx)
DEF_MEMx(stwux)
DEF_MEMx(stwbrx)
DEF_MEMx(sthx)
DEF_MEMx(sthbrx)
DEF_MEMx(stdx)
DEF_MEMx(stdcx)
DEF_MEMx(stdux)
DEF_MEMx(stwcx)
#undef DEF_MEMx

BufferOffset Assembler::as_isel(Register rt, Register ra, Register rb, uint32_t bc)
{
    spew("isel\t%3s,%3s,%3s,%d", rt.name(), ra.name(), rb.name(), bc);
    return writeInst(PPC_isel | rt.code() << 21 | ra.code() << 16 | rb.code() << 11 | bc << 6);
}

// FPR operations and load-stores.
BufferOffset Assembler::as_fcmpo(CRegisterID cr, FloatRegister ra, FloatRegister rb)
{
    spew("fcmpo\t%d,%3s,%3s", cr, ra.name(), rb.name());
    return writeInst(PPC_fcmpo | cr << 23 | ra.code() << 16 | rb.code() << 11);
}

BufferOffset Assembler::as_fcmpo(FloatRegister ra, FloatRegister rb)
{
    return as_fcmpo(cr0, ra, rb);
}

BufferOffset Assembler::as_fcmpu(CRegisterID cr, FloatRegister ra, FloatRegister rb)
{
    spew("fcmpu\t%d,%3s,%3s", cr, ra.name(), rb.name());
    return writeInst(PPC_fcmpu | cr << 23 | ra.code() << 16 | rb.code() << 11);
}

BufferOffset Assembler::as_fcmpu(FloatRegister ra, FloatRegister rb)
{
    return as_fcmpu(cr0, ra, rb);
}

#define DEF_FPUAC(op) DEF_AFORM_C(op) DEF_AFORM_C_RC(op)
DEF_FPUAC(fmul)
DEF_FPUAC(fmuls)
#undef DEF_FPUAC

#define DEF_FPUAB(op) DEF_AFORM_B(op) DEF_AFORM_B_RC(op)
DEF_FPUAB(fadd)
DEF_FPUAB(fdiv)
DEF_FPUAB(fsub)
DEF_FPUAB(fadds)
DEF_FPUAB(fdivs)
DEF_FPUAB(fsubs)
DEF_FPUAB(fcpsgn)
#undef DEF_FPUAB

#define DEF_FPUDS(op) DEF_XFORM2_F(op) DEF_XFORM2_F_RC(op)
DEF_FPUDS(fabs)
DEF_FPUDS(fneg)
DEF_FPUDS(fmr)
DEF_FPUDS(fcfid)
DEF_FPUDS(fcfidu)
DEF_FPUDS(fctid)
//DEF_FPUDS(fctidz)
DEF_FPUDS(fctidu)
//DEF_FPUDS(fctiduz)
DEF_FPUDS(fctiw)
DEF_FPUDS(fctiwz)
DEF_FPUDS(fctiwu)
DEF_FPUDS(fctiwuz)
DEF_FPUDS(frim)
DEF_FPUDS(frin)
DEF_FPUDS(frip)
DEF_FPUDS(friz)
DEF_FPUDS(frsp)
DEF_FPUDS(frsqrte)
DEF_FPUDS(fsqrt)
DEF_FPUDS(fsqrts)
#undef DEF_FPUDS

// In Ion, the semantics for this macro are now corrected compared to JM/PPCBC. 
// (See OPPCC p.432, etc.)
#define DEF_FPUACB(op) DEF_AFORM(op) DEF_AFORM_RC(op)
DEF_FPUACB(fmadd)
DEF_FPUACB(fnmsub)
DEF_FPUACB(fsel)
#undef DEF_FPUACB

// XXX: This displays offsets wrong in spew (we want f,d(r))
#define DEF_FMEMd(op) DEF_DFORM_F(op)
DEF_FMEMd(lfd)
DEF_FMEMd(lfs)
DEF_FMEMd(stfd)
DEF_FMEMd(stfs)
DEF_FMEMd(stfdu)
DEF_FMEMd(stfsu)
#undef DEF_FMEMd

#define DEF_FMEMx(op) \
    BufferOffset Assembler::as_##op(FloatRegister rd, Register ra, Register rb) { \
        spew(#op "\t%3s,%3s,%3s", rd.name(), ra.name(), rb.name()); \
        return writeInst(XForm(PPC_##op, rd, ra, rb, false)); }
DEF_FMEMx(lfdx)
DEF_FMEMx(lfsx)
DEF_FMEMx(lfiwax)
DEF_FMEMx(stfiwx)
DEF_FMEMx(stfdx)
DEF_FMEMx(stfsx)
#undef DEF_FMEMx

BufferOffset Assembler::as_mtfsb0(uint8_t bt)
{
    spew("mtfsb0\t%d", bt);
    return writeInst(PPC_mtfsb0 | (uint32_t)bt << 21);
}

BufferOffset Assembler::as_mtfsb1(uint8_t bt)
{
    spew("mtfsb1\t%d", bt);
    return writeInst(PPC_mtfsb1 | (uint32_t)bt << 21);
}

BufferOffset Assembler::as_mtfsfi(uint8_t fi, uint8_t imm)
{
    spew("mtfsfi\t%d,%d", fi, imm);
    return writeInst(PPC_mtfsfi | fi << 23 | imm << 12);
}

BufferOffset Assembler::as_mcrf(CRegisterID bt, CRegisterID bs)
{
    spew("mcrf\t%d,%d", bt, bs);
    return writeInst(PPC_mcrf | (uint32_t)bt << 23 | (uint32_t)bs << 18);
}

BufferOffset Assembler::as_mcrfs(CRegisterID bf, uint8_t bfa)
{
    spew("mcrfs\t%d,%d", bf, bfa);
    return writeInst(PPC_mcrfs | (uint32_t)bf << 23 | (uint32_t)bfa << 18);
}

// VSX
// No RC forms for these (least significant bit sets vector or FPR).
    BufferOffset Assembler::as_mfvsrd(Register ra, FloatRegister xs) {
        spew("mfvsrd\t%3s,%3s", ra.name(), xs.name());
        return writeInst(XForm(PPC_mfvsrd, xs, ra, r0, false));
    }
    BufferOffset Assembler::as_mtvsrd(FloatRegister xt, Register ra) {
        spew("mtvsrd\t%3s,%3s", xt.name(), ra.name());
        // Yes, same operand order (see PowerISA v3.1 page 121)
        return writeInst(XForm(PPC_mtvsrd, xt, ra, r0, false));
    }

// Conveniences and generally accepted alternate mnemonics.
// XXX: change these to xs_
BufferOffset Assembler::xs_trap()
{
    spew("trap");
    return writeInst(PPC_trap);
}
// trap with metadata encoded as register
BufferOffset Assembler::xs_trap_tagged(uint8_t tag) {
    // FreeBSD and others may use r1 in their trap word, so don't allow bit 0 or > 15.
    spew("trap ; MARK %d", tag);
    MOZ_ASSERT(!(tag & 1) && (tag < 16));
    return writeInst(PPC_trap | (tag << 21) | (tag << 16));
}

BufferOffset Assembler::x_mr(Register rd, Register ra)
{
    return as_or(rd, ra, ra);
}

BufferOffset Assembler::x_beq(CRegisterID cr, int16_t off, LikelyBit lkb, LinkBit lb)
{
    return as_bc(off, Equal, cr, lkb, lb);
}

BufferOffset Assembler::x_bne(CRegisterID cr, int16_t off, LikelyBit lkb, LinkBit lb)
{
    return as_bc(off, NotEqual, cr, lkb, lb);
}

BufferOffset Assembler::x_bdnz(int16_t off, LikelyBit lkb, LinkBit lb)
{
    spew("bdnz %d\n", off);
    return writeInst(PPC_bc | 0x10 << 21 | off | lkb << 21 | lb);
}

BufferOffset Assembler::xs_mtctr(Register ra)
{
    return as_mtspr(ctr, ra);
}

BufferOffset Assembler::xs_mtlr(Register ra)
{
    return as_mtspr(lr_spr, ra);
}

BufferOffset Assembler::xs_mflr(Register rd)
{
    return as_mfspr(rd, lr_spr);
}

BufferOffset Assembler::xs_mtcr(Register rs)
{
    return as_mtcrf(0xff, rs);
}

BufferOffset Assembler::xs_mfxer(Register ra)
{
    return as_mfspr(ra, xer);
}

BufferOffset Assembler::xs_mtxer(Register ra)
{
    return as_mtspr(xer, ra);
}

BufferOffset Assembler::x_bit_value(Register rd, Register rs, unsigned bit)
{
    return as_rlwinm(rd, rs, bit + 1, 31, 31);
}

BufferOffset Assembler::x_slwi(Register rd, Register rs, int n)
{
    return as_rlwinm(rd, rs, n, 0, 31 - n);
}

BufferOffset Assembler::x_sldi(Register rd, Register rs, int n)
{
    return as_rldicr(rd, rs, n, 63 - n);
}

BufferOffset Assembler::x_srwi(Register rd, Register rs, int n)
{
    return as_rlwinm(rd, rs, 32 - n, n, 31);
}

BufferOffset Assembler::x_srdi(Register rd, Register rs, int n)
{
    return as_rldicl(rd, rs, 64 - n, n);
}

BufferOffset Assembler::x_subi(Register rd, Register ra, int16_t im)
{
    return as_addi(rd, ra, -im);
}

BufferOffset Assembler::xs_li(Register rd, int16_t im)
{
    return as_addi(rd, r0, im, true /* actually_li */);
}

BufferOffset Assembler::xs_lis(Register rd, int16_t im)
{
    return as_addis(rd, r0, im, true /* actually_lis */);
}

BufferOffset Assembler::x_not(Register rd, Register ra)
{
    return as_nor(rd, ra, ra);
}

// Traps
BufferOffset Assembler::as_tw(uint8_t to, Register ra, Register rb)
{
    return writeInst(PPC_tw | (uint32_t)to << 21 | ra.code() << 16 |
            rb.code() << 11);
}

BufferOffset Assembler::as_twi(uint8_t to, Register ra, int16_t si)
{
    return writeInst(PPC_twi | (uint32_t)to << 21 | ra.code() << 16 | si);
}
