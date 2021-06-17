/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef jit_ppc_Architecture_ppc_h
#define jit_ppc_Architecture_ppc_h

#include <limits.h>
#include <stdint.h>

#include "mozilla/Assertions.h"
#include "mozilla/MathAlgorithms.h"

#include "jit/shared/Architecture-shared.h"

#include "js/Utility.h"

namespace js {
namespace jit {

// Not used on PPC.
static const uint32_t ShadowStackSpace = 0;

// Size of each bailout table entry.
// For PowerPC this is a single bl.
static const uint32_t BAILOUT_TABLE_ENTRY_SIZE = sizeof(void *);

// Range of an immediate jump (26 bit jumps
static constexpr uint32_t JumpImmediateRange = 32 * 1024 * 1024;

// GPRs.
class Registers
{
  public:
    enum RegisterID {
        r0 = 0,
        tempRegister = r0,
        r1,
        sp = r1,
        stackPointerRegister = r1,
        r2,
        r3,
        r4,
        r5,
        r6,
        r7,
        r8,
        r9,
        r10,
        r11,
        r12,
        addressTempRegister = r12,
        r13,
        r14,
        r15,
        r16,
        r17,
        r18,
        r19,
        r20,
        r21,
        r22,
        r23,
        r24,
        r25,
        r26,
        r27,
        r28,
        r29,
        r30,
        r31,
        invalid_reg
    };

    typedef uint8_t Code;
    typedef uint32_t Encoding;
    typedef uint64_t SetType;
    
    // Content spilled during bailouts.
    union RegisterContent {
    	uintptr_t r;
    };

    static const char *GetName(Code code) {
        static const char *Names[] = {
             "r0",  "sp",  "toc", "r3",  "r4",  "r5",  "r6",  "r7",
             "r8",  "r9",  "r10", "r11", "r12", "r13", "r14", "r15",
             "r16", "r17", "r18", "r19", "r20", "r21", "r22", "r23",
             "r24", "r25", "r26", "r27", "r28", "r29", "r30", "r31"};
        return Names[code];
    }
    static const char *GetName(uint32_t i) {
        MOZ_ASSERT(i < Total);
        return GetName(Code(i));
    }

    static Code FromName(const char *name);

    static const Encoding StackPointer = sp;
    static const Encoding Invalid = invalid_reg;

    // If we decide to implement any virtual registers or subsume SPRs
    // into the GPR space, we have the headroom for it.
    static const uint32_t Total = 32;
    static const uint32_t Allocatable = 24;

    static const SetType AllMask = 0x00000000ffffffff;
    static const SetType ArgRegMask =
        (1 << Registers::r3) |
        (1 << Registers::r4) |
        (1 << Registers::r5) |
        (1 << Registers::r6) |
        (1 << Registers::r7) |
        (1 << Registers::r8) |
        (1 << Registers::r9) |
        (1 << Registers::r10);

    // Don't bother saving r0, r1, r11 or r12.
    static const SetType VolatileMask = ArgRegMask;

    // We use this constant to save registers when entering functions.
    static const SetType NonVolatileMask = (
    	(1 << Registers::r2)  |
        (1 << Registers::r13) |
        (1 << Registers::r14) |
        (1 << Registers::r15) |
        (1 << Registers::r16) |
        (1 << Registers::r17) |
        (1 << Registers::r18) |
        (1 << Registers::r19) |
        (1 << Registers::r20) |
        (1 << Registers::r21) |
        (1 << Registers::r22) |
        (1 << Registers::r23) |
        (1 << Registers::r24) |
        (1 << Registers::r25) |
        (1 << Registers::r26) |
        (1 << Registers::r27) |
        (1 << Registers::r28) |
        (1 << Registers::r29) |
        (1 << Registers::r30) |
        (1 << Registers::r31)
    // Watch out for sign extension!
        ) & AllMask;

    // Also uses r11.
    static const SetType WrapperMask = VolatileMask;

    static const SetType NonAllocatableMask =
        // Used by assembler.
        (1 << Registers::r0)  |
        (1 << Registers::sp)  |
        (1 << Registers::r2)  |
        // Temp registers.
        (1 << Registers::r11) |
        (1 << Registers::r12) |
        // r13 is the pointer for TLS in ELF v2.
        (1 << Registers::r13) |
        // Non-volatile work registers.
        (1 << Registers::r16) |
        // r17 is the InterpreterPCReg and must be allocatable.
        (1 << Registers::r18);
        // Despite its use as a rectifier, r19 must be allocatable (see
        // ICCallScriptedCompiler::generateStubCode).

    // Registers that can be allocated without being saved, generally.
    static const SetType TempMask = VolatileMask & ~NonAllocatableMask;

    // Registers returned from a JS -> JS call.
    static const SetType JSCallMask =
        (1 << Registers::r5);

    // Registers returned from a JS -> C call.
    static const SetType CallMask =
        (1 << Registers::r3);
 
    static const SetType AllocatableMask = (
    	// Be explicit
        (1 << Registers::r3)  |
        (1 << Registers::r4)  |
        (1 << Registers::r5)  |
        (1 << Registers::r6)  |
        (1 << Registers::r7)  |
        //(1 << Registers::r8)  |
        (1 << Registers::r9)  |
        (1 << Registers::r10) |
        (1 << Registers::r14) |
        (1 << Registers::r15) |
        (1 << Registers::r17) |
        (1 << Registers::r19) |
        (1 << Registers::r20) |
        (1 << Registers::r21) |
        (1 << Registers::r22) |
        (1 << Registers::r23) |
        (1 << Registers::r24) |
        (1 << Registers::r25) |
        (1 << Registers::r26) |
        (1 << Registers::r27) |
        (1 << Registers::r28) |
        (1 << Registers::r29) |
        (1 << Registers::r30) |
        (1 << Registers::r31)
        // Watch out for sign extension!
        ) & AllMask;

    static uint32_t SetSize(SetType x) {
        static_assert(sizeof(SetType) == 8, "SetType must be 64 bits");
        return mozilla::CountPopulation64(x);
    }
    static uint32_t FirstBit(SetType x) {
        return mozilla::CountTrailingZeroes64(x);
    }
    static uint32_t LastBit(SetType x) {
        return 63 - mozilla::CountLeadingZeroes64(x);
    }
};

// Smallest integer type that can hold a register bitmask.
typedef uint64_t PackedRegisterMask;

// FPRs.
// PowerPC FPRs can be both double and single precision, like MIPS. We tell
// Ion there are 64 FPRs, but each is an aliased pair.
class FloatRegisters
{
  public:
    enum FPRegisterID {
        f0 = 0,
        f1,
        f2,
        f3,
        f4,
        f5,
        f6,
        f7,
        f8,
        f9,
        f10,
        f11,
        f12,
        f13,
        f14,
        f15,
        f16,
        f17,
        f18,
        f19,
        f20,
        f21,
        f22,
        f23,
        f24,
        f25,
        f26,
        f27,
        f28,
        f29,
        f30,
        f31,
        invalid_freg
    };
    typedef uint32_t Code;
    typedef uint64_t SetType;
    typedef FPRegisterID Encoding;
    enum ContentType { Single, Double, NumTypes };

  // Content spilled during bailouts.
  union RegisterContent {
    float s;
    double d;
  };

  // Only the double registers (i.e., the underlying FPRs) need actually be saved, not aliases.
  static const uint32_t TotalPhys = 32;
  static const uint32_t Total = TotalPhys * NumTypes;
  //static const uint32_t TotalWithSimd = TotalPhys;
  static const uint32_t Allocatable = 58;  // Don't use f0, f2 or f3.

  static_assert(sizeof(SetType) * 8 >= Total,
                "SetType should be large enough to enumerate all registers.");

  // Magic values which are used to duplicate a mask of physical registers for
  // a specific type of register. A multiplication is used to copy and shift
  // the bits of the physical register mask.
  static const SetType SpreadSingle = SetType(1)
                                      << (uint32_t(Single) * TotalPhys);
  static const SetType SpreadDouble = SetType(1)
                                      << (uint32_t(Double) * TotalPhys);
  static const SetType SpreadScalar = SpreadSingle | SpreadDouble;
  static const SetType SpreadVector = 0;
  static const SetType Spread = SpreadScalar | SpreadVector;

  static const SetType AllPhysMask = ((SetType(1) << TotalPhys) - 1);
  static const SetType AllMask = AllPhysMask * Spread;
  static const SetType AllSingleMask = AllPhysMask * SpreadSingle;
  static const SetType AllDoubleMask = AllPhysMask * SpreadDouble;

   static const char *GetName(Encoding code) {
        static const char * const Names[] = { "f0", "f1", "f2", "f3",  "f4", "f5",  "f6", "f7",
                                              "f8", "f9",  "f10", "f11", "f12", "f13",
                                              "f14", "f15", "f16", "f17", "f18", "f19",
                                              "f20", "f21", "f22", "f23", "f24", "f25",
                                              "f26", "f27", "f28", "f29", "f30", "f31"};
        MOZ_ASSERT(code < TotalPhys);
        return Names[code];
    }
    static const char *GetName(uint32_t i) {
        MOZ_ASSERT(i < Total);
        return GetName(Code(i % TotalPhys));
    }

    static Encoding FromName(const char *name);

    static const Encoding Invalid = invalid_freg;

    static const SetType VolatileMask = (
        // Since we never let it allocate f0, don't let it push it either.
        (1 << FloatRegisters::f1)  |
        (1 << FloatRegisters::f2)  |
        (1 << FloatRegisters::f3)  |
        (1 << FloatRegisters::f4)  |
        (1 << FloatRegisters::f5)  |
        (1 << FloatRegisters::f6)  |
        (1 << FloatRegisters::f7)  |
        (1 << FloatRegisters::f8)  |
        (1 << FloatRegisters::f9)  |
        (1 << FloatRegisters::f10) |
        (1 << FloatRegisters::f11) |
        (1 << FloatRegisters::f12) |
        (1 << FloatRegisters::f13)) * SpreadScalar | AllPhysMask * SpreadVector;

    static const SetType NonVolatileMask = (
        (1 << FloatRegisters::f14) |
        (1 << FloatRegisters::f15) |
        (1 << FloatRegisters::f16) |
        (1 << FloatRegisters::f17) |
        (1 << FloatRegisters::f18) |
        (1 << FloatRegisters::f19) |
        (1 << FloatRegisters::f20) |
        (1 << FloatRegisters::f21) |
        (1 << FloatRegisters::f22) |
        (1 << FloatRegisters::f23) |
        (1 << FloatRegisters::f24) |
        (1 << FloatRegisters::f25) |
        (1 << FloatRegisters::f26) |
        (1 << FloatRegisters::f27) |
        (1 << FloatRegisters::f28) |
        (1 << FloatRegisters::f29) |
        (1 << FloatRegisters::f30) |
        (1 << FloatRegisters::f31)) * SpreadScalar | AllPhysMask * SpreadVector;

	//static const SetType VolatileDoubleMask = VolatileMask;
	//static const SetType NonVolatileDoubleMask = NonVolatileMask;

    static const SetType WrapperMask = VolatileMask;

    // f1 must be allocatable because Ion may expect to optimize with it.
    static const SetType NonAllocatableMask = (
        (1 << FloatRegisters::f0) |
        (1 << FloatRegisters::f2) |
        (1 << FloatRegisters::f3)) * Spread;
	//static const SetType NonAllocatableDoubleMask = NonAllocatableMask;

    // Registers that can be allocated without being saved, generally.
    //static const SetType TempMask = VolatileMask & ~NonAllocatableMask;
	
    static const SetType AllocatableMask = (
        // Be explicit
        (1 << FloatRegisters::f1)  |
        (1 << FloatRegisters::f4)  |
        (1 << FloatRegisters::f5)  |
        (1 << FloatRegisters::f6)  |
        (1 << FloatRegisters::f7)  |
        (1 << FloatRegisters::f8)  |
        (1 << FloatRegisters::f9)  |
        (1 << FloatRegisters::f10) |
        (1 << FloatRegisters::f11) |
        (1 << FloatRegisters::f12) |
        (1 << FloatRegisters::f13) |
        (1 << FloatRegisters::f14) |
        (1 << FloatRegisters::f15) |
        (1 << FloatRegisters::f16) |
        (1 << FloatRegisters::f17) |
        (1 << FloatRegisters::f18) |
        (1 << FloatRegisters::f19) |
        (1 << FloatRegisters::f20) |
        (1 << FloatRegisters::f21) |
        (1 << FloatRegisters::f22) |
        (1 << FloatRegisters::f23) |
        (1 << FloatRegisters::f24) |
        (1 << FloatRegisters::f25) |
        (1 << FloatRegisters::f26) |
        (1 << FloatRegisters::f27) |
        (1 << FloatRegisters::f28) |
        (1 << FloatRegisters::f29) |
        (1 << FloatRegisters::f30) |
        (1 << FloatRegisters::f31)) * Spread;
    //static const SetType AllocatableDoubleMask = AllocatableMask;
};

template <typename T>
class TypedRegisterSet;

// Shamelessly rip off MIPS
class FloatRegister
{
 public:
  typedef FloatRegisters Codes;
  typedef size_t Code;
  typedef Codes::Encoding Encoding;
  typedef Codes::ContentType ContentType;
  typedef FloatRegisters::SetType SetType;

  Encoding reg_ : 6;

 private:
  ContentType kind_ : 3;

 public:
  constexpr FloatRegister(uint32_t r, ContentType kind = Codes::Double)
      : reg_(Encoding(r)), kind_(kind) {}
  constexpr FloatRegister()
      : reg_(Encoding(FloatRegisters::invalid_freg)), kind_(Codes::Double) {}

  static uint32_t SetSize(SetType x) {
    // Count the number of non-aliased registers.
    x |= x >> Codes::TotalPhys;
    x &= Codes::AllPhysMask;
    return mozilla::CountPopulation64(x);
  }
  static uint32_t FirstBit(SetType x) {
    static_assert(sizeof(SetType) == 8, "SetType must be 64 bits");
    return mozilla::CountTrailingZeroes64(x);
  }
  static uint32_t LastBit(SetType x) {
    static_assert(sizeof(SetType) == 8, "SetType must be 64 bits");
    return 63 - mozilla::CountLeadingZeroes64(x);
  }

  bool operator==(const FloatRegister& other) const {
    MOZ_ASSERT(!isInvalid());
    MOZ_ASSERT(!other.isInvalid());
    return kind_ == other.kind_ && reg_ == other.reg_;
  }
  bool equiv(const FloatRegister& other) const { return other.kind_ == kind_; }
  size_t size() const {
    return (kind_ == Codes::Double) ? sizeof(double) : sizeof(float);
  }
  // Always push doubles to maintain 8-byte stack alignment.
  size_t pushSize() const { return sizeof(double); }
  bool isInvalid() const { return reg_ == FloatRegisters::invalid_freg; }

  bool isSingle() const { return kind_ == Codes::Single; }
  bool isDouble() const { return kind_ == Codes::Double; }
  bool isSimd128() const { return false; }

  FloatRegister singleOverlay() const;
  FloatRegister doubleOverlay() const;

  FloatRegister asSingle() const { return singleOverlay(); }
  FloatRegister asDouble() const { return doubleOverlay(); }
  FloatRegister asSimd128() const { MOZ_CRASH("NYI"); }

  Code code() const {
    MOZ_ASSERT(!isInvalid());
    return Code(reg_ | (kind_ << 5));
  }
  Encoding encoding() const {
    MOZ_ASSERT(!isInvalid());
    MOZ_ASSERT(uint32_t(reg_) < Codes::TotalPhys);
    return reg_;
  }
  uint32_t id() const { return reg_; }
  static FloatRegister FromCode(uint32_t i) {
    uint32_t code = i & 0x1f;
    uint32_t kind = i >> 5;
    return FloatRegister(Code(code), ContentType(kind));
  }

  bool volatile_() const {
    return !!((1 << reg_) & FloatRegisters::VolatileMask);
  }
  const char* name() const { return FloatRegisters::GetName(reg_); }
  bool operator!=(const FloatRegister& other) const {
    return kind_ != other.kind_ || reg_ != other.reg_;
  }
  bool aliases(const FloatRegister& other) { return reg_ == other.reg_; }
  uint32_t numAliased() const { return 2; }
  FloatRegister aliased(uint32_t aliasIdx) {
    if (aliasIdx == 0) {
      return *this;
    }
    MOZ_ASSERT(aliasIdx == 1);
    if (isDouble()) {
      return singleOverlay();
    }
    return doubleOverlay();
  }
  uint32_t numAlignedAliased() const { return 2; }
  FloatRegister alignedAliased(uint32_t aliasIdx) {
    MOZ_ASSERT(isDouble());
    if (aliasIdx == 0) {
      return *this;
    }
    MOZ_ASSERT(aliasIdx == 1);
    return singleOverlay();
  }

  SetType alignedOrDominatedAliasedSet() const { return Codes::Spread << reg_; }

  static constexpr RegTypeName DefaultType = RegTypeName::Float64;

  template <RegTypeName = DefaultType>
  static SetType LiveAsIndexableSet(SetType s) {
    return SetType(0);
  }
  template <RegTypeName Name = DefaultType>
  static SetType AllocatableAsIndexableSet(SetType s) {
    static_assert(Name != RegTypeName::Any, "Allocatable set are not iterable");
    return LiveAsIndexableSet<Name>(s);
  }

  static Code FromName(const char* name) {
    return FloatRegisters::FromName(name);
  }
  static TypedRegisterSet<FloatRegister> ReduceSetForPush(
      const TypedRegisterSet<FloatRegister>& s);
  static uint32_t GetPushSizeInBytes(const TypedRegisterSet<FloatRegister>& s);
  uint32_t getRegisterDumpOffsetInBytes();
};

template <>
inline FloatRegister::SetType
FloatRegister::LiveAsIndexableSet<RegTypeName::Float32>(SetType set) {
  return set & FloatRegisters::AllSingleMask;
}

template <>
inline FloatRegister::SetType
FloatRegister::LiveAsIndexableSet<RegTypeName::Float64>(SetType set) {
  return set & FloatRegisters::AllDoubleMask;
}

template <>
inline FloatRegister::SetType
FloatRegister::LiveAsIndexableSet<RegTypeName::Any>(SetType set) {
  return set;
}

// XXX: This needs to be rolled into somewhere else

// SPRs (PPC backend specific).
// These have no peer in lesser chips. That is because PPC has no peer in
// lesser chips. These don't count against the register cap because the
// allocator is unaware of them. In fact, we don't treat these as regular
// registers at all (hey, they're Special Purpose anyway).
#if(0)
class SPRs // Class definition not yet supported.
{
  public:
#endif
    enum SPRegisterID {
      xer = 1,
      lr_spr = 8, /* I'm in the REAL witness protection program. */
      ctr = 9,
      vrsave = 256, /* for future SIMD JS */
      invalid_spreg
    };
#if(0)
    static const char *getSPRName(SPRegisterID code) {
#define XXX "INVALID"
        static const char *N_vrsave = "vrsave";
        static const char *N_bogus = XXX;
        static const char *Names[] = {
            XXX, "xer", XXX, XXX, XXX, XXX, XXX, XXX,
            "lr","ctr"
        };
#undef XXX
        return 
               (code == vrsave) ? N_vrsave :
               (code >  ctr)    ? N_bogus :
               Names[code];
    }
};
#endif

// CRs (PPC backend specific).
// We have eight condition registers, each for how unconditionally wonderful
// PowerPC is, and sometimes for storing condition results.
// Assume CR0 as default.
#if(0)
class CRs // Class definition not yet supported.
{
  public:
#endif
    enum CRegisterID {
      cr0 = 0,
      cr1,
/*
Suppress CR2-CR4 because these are non-volatile.
      cr2,
      cr3,
      cr4,
*/
      cr5 = 5,
      cr6,
      cr7,
      invalid_creg
    };
#if(0)
    static const char *getCRName(CRegisterID code) {
        static const char *Names[] = {
            "cr0",  "cr1",  "cr2",  "cr3",  "cr4",  "cr5",  "cr6",  "cr7"
        };
        return Names[code];
    }
#endif

inline uint32_t GetPPC64Flags() { return 0; /* XXX? */ }
inline bool hasUnaliasedDouble() { return false; }
inline bool hasMultiAlias() { return false; }

} // namespace jit
} // namespace js

#endif /* jit_ppc_Architecture_ppc_h */
