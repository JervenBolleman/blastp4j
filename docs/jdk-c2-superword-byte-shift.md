# C2 SuperWord: `short[i] = (short) (byteArray[i] << k)` is vectorized through a byte and loses data

A loop that reads a `byte[]`, shifts left, and stores into a `short[]` or `char[]` is auto-vectorized into
code that truncates each result to 8 bits and sign-extends it back to 16, so every value that does not fit
a byte comes out wrong. Scalar execution, `-XX:-UseSuperWord` and `-XX:TieredStopAtLevel=1` are all correct.

## Affected

| JDK | result |
|---|---|
| 21.0.13.1 (Corretto) | OK |
| 25.0.1 (OpenJDK) | OK |
| 25.0.5 (OpenJDK), 25.0.5.1 (Corretto) | OK |
| **27+35-2325 (OpenJDK, 2026-09-15)** | **wrong** |
| **27.0.1.0.3 (OpenJDK, 2026-10-24)** | **wrong** |

No JDK 26 build was available to narrow the window further.

Platform: Intel Core i7-6700HQ (Skylake, AVX2, no AVX-512), Fedora 45, Linux 7.3.4-300.fc45.x86_64.
Reproduces with no flags at all.

## Reproducer

`ShiftFill.java` (also in this directory, with `ShiftFillChar.java` and `ShiftFillInt.java`):

```java
public class ShiftFill {

	static short[] fill(byte[] in) {
		short[] out = new short[in.length];
		for (int i = 0; i < in.length; i++) {
			out[i] = (short) (in[i] << 5);
		}
		return out;
	}

	public static void main(String[] args) {
		final int n = 142, rounds = 200_000;
		byte[] in = new byte[n];
		for (int i = 0; i < n; i++) {
			in[i] = (byte) (i % 25);
		}
		short[] want = new short[n];                 // computed before anything is hot
		for (int i = 0; i < n; i++) {
			want[i] = (short) (in[i] << 5);
		}
		int bad = 0;
		String first = null;
		for (int r = 0; r < rounds; r++) {
			short[] got = fill(in);
			for (int i = 0; i < n; i++) {
				if (got[i] != want[i]) {
					bad++;
					if (first == null) {
						first = "round " + r + ", index " + i + ": in=" + in[i]
								+ " want=" + want[i] + " got=" + got[i];
					}
					break;
				}
			}
		}
		System.out.println(bad == 0 ? "OK" : "WRONG in " + bad + " of " + rounds + " rounds; first " + first);
	}
}
```

```
$ javac ShiftFill.java && java ShiftFill
WRONG in 185581 of 200000 rounds, n=142; first round 14419, index 6: in=6 want=192 got=-64

$ java -XX:-UseSuperWord ShiftFill
OK (200000 rounds, n=142)

$ java -XX:TieredStopAtLevel=1 ShiftFill
OK (200000 rounds, n=142)
```

`6 << 5` is 192; the vectorized loop returns −64, which is 192 truncated to a byte (`0xC0`) and
sign-extended back to 16 bits (`0xFFC0`). Inputs 0..3 give 0, 32, 64 and 96, which fit a byte and survive;
every input from 4 up is corrupted. The first wrong index varies between runs because the scalar prologue
covers a different number of leading elements.

`ShiftFillChar.java` is the same with `char[]`: `in=4`, want 128, got 65408 (`0xFF80`) — the same byte
round-trip, so the bug is in the 16-bit narrowing store rather than in `short` specifically.
`ShiftFillInt.java` stores into an `int[]`, has no narrowing store, and is correct.

## Generated code

`-XX:CompileCommand=print,ShiftFill::fill`, disassembled with objdump (no hsdis on this machine):

```asm
vmovq      0xc(%rdx,%rsi,1),%xmm2    ; load 8 bytes of in[]
vpmovsxbw  %xmm2,%xmm6               ; byte -> word, sign-extended        correct
vpsllw     %xmm0,%xmm6,%xmm6         ; shift the 8 words left by 5        correct
vmovdqu    -0x523473(%rip),%xmm5     ; a per-lane 0x00ff mask
vpand      %xmm6,%xmm5,%xmm5         ; keep the low 8 bits of each word   <-- loses bit 7 upward
vpackuswb  %xmm5,%xmm5,%xmm5         ; pack the words down to bytes
vpmovsxbw  %xmm5,%xmm2               ; sign-extend the bytes back to words
vmovdqu    %xmm2,0xc(%r13,%r10,2)    ; store 8 shorts
```

The value is computed correctly by `vpsllw` and then round-tripped through a byte by the
`vpand` / `vpackuswb` / `vpmovsxbw` triple. The source casts to `short`, so the generated code is
computing `(short) (byte) (in[i] << 5)`.

The narrowing looks like it is carried over from the element type of the *input* array: the shift widens
the value range beyond a byte, and that seems not to be accounted for when the pack is chosen.

## Workarounds

- `-XX:-UseSuperWord` (global, and costly elsewhere).
- Avoid the narrowing store: keep the destination an `int[]`, or do the shift at the point of use.

## How it was found

In blastp4j (https://github.com/JervenBolleman/blastp4j), narrowing a `rowOfQuery` table from `int[]` to
`short[]` made a protein search silently miss alignments. Through the application the failure looked
flaky — about 5% of test runs — because it depended on whether the constructor holding the loop had been
C2-compiled yet; instances built before compilation had correct tables and instances built after did not.
Driving the method hot directly, as the reproducer above does, makes it deterministic.

## Where it probably comes from

`VLoopTypes::compute_vector_element_type` in `src/hotspot/share/opto/superword.cpp` propagates narrowed
integer types backwards through the loop body, gated by `can_subword_truncate`. That gate was added by
[openjdk/jdk#25440](https://github.com/openjdk/jdk/pull/25440),
[JDK-8350177](https://bugs.openjdk.org/browse/JDK-8350177), integrated 2025-07-13 — after the JDK 25 fork,
which fits the regression window above. Its own description says it fixes "cases in SuperWord when
compiling subword types where vectorized code would be given a narrower type than expected, leading to
miscompilation due to truncation", and that full correctness waits on
[JDK-8342095](https://bugs.openjdk.org/browse/JDK-8342095) adding automatic cast insertion. A steady
stream of follow-ups has landed since (JDK-8362171, 8362972, 8362979, 8365570, 8369881, 8380988, 8394230,
8394456), all of the form "unexpected node in SuperWord truncation".

The relevant lines carry an assumption this case breaks:

```cpp
if ((VectorNode::is_shift_opcode(op) || op == Op_AbsI) && load->is_Load() && ...) {
    vt = velt_type(load);        // for shifts, take the element type from the input load
} else if (op != Op_LShiftI) {
    // Widen type to the node type to avoid the creation of vector nodes. Note
    // that left shifts work regardless of the signedness.
    vt = container_type(in);
}
```

`Op_LShiftI` is exempted from widening because a left shift is insensitive to the signedness of its input.
That is true, but it is not the property needed here: the result of this shift is *stored wider than its
input was loaded*, so an element type derived from the `byte[]` load discards exactly the bits the shift
moved up. `can_subword_truncate` also returns `true` for any `is_Load()` without qualification.

This is a reading of the source rather than a diagnosis from a debug build, so it is offered as a starting
point, not a conclusion.
