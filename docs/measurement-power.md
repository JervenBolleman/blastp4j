# How many rounds a verdict needs

Every experiment in this repository ends in a verdict, and most of them have ended in "flat".
This note works out what that word is allowed to mean, from the measurements the runners have
already made, and how many rounds it takes to earn it.

## "Flat" is a claim, not the absence of one

A difference that does not reach significance is not evidence of no difference. The two are
distinguishable only by the width of the interval around the estimate: a result can fail to be
significant because the effect is small, or because the measurement is too noisy to see it, and
nothing in the p-value tells the two apart.

Claiming flatness is therefore an equivalence claim, and it needs a margin. The right form is
two one-sided tests: fix a margin `d` beforehand, and declare the change flat when the whole
99% confidence interval of the difference lies inside `(-d, +d)`. Without a stated `d` the word
has no content, because no experiment can ever show an effect is exactly zero.

The margin is a judgement about what matters, not a statistic. For this project 1.5% is a
defensible choice: the spread of true effects across the experiments measured so far is about
1.75%, so a 1.5% margin rules out everything worth having while staying reachable.

## What the noise actually is

From the 88 cell-measurements of the human-proteome workload (`UP000005640_9606`), each cell
being one architecture, one program and three rounds of the candidate, the reference and BLAST+:

| | round-to-round CV |
|---|---|
| BLAST+, Neoverse-N2 full | 0.34% |
| blastp4j, Neoverse-N2 full | 1.22% candidate, 1.61% reference |
| BLAST+, Xeon 6975P-C full | 0.42% |
| blastp4j, Xeon 6975P-C full | 1.05% candidate, 1.54% reference |

The machine is quiet. BLAST+ is the same binary on the same host doing the same work, and it
repeats to a third of a percent. The JVM on that host repeats four times worse. The noise this
project fights is the JVM's own, not the runner's.

Three further properties, all measured rather than assumed:

- **The run-level effect cancels.** BLAST+ varies between runs by 0.47% beyond what its rounds
  explain, so runs do differ. But the candidate and the reference are measured inside the same
  run, and that offset cancels in their ratio. Predicted from round noise alone, the paired
  statistic should have sd 1.17%; measured over ten cells where the two commits compile
  identical machine code, it is **1.13%**. There is no residual run-level term to fight.
- **Rounds are exchangeable.** Round 1, 2 and 3 average -0.01%, +0.02% and -0.01% of their
  cell's mean. There is no cold first round to discard.
- **The noise is symmetric.** The mean exceeds the median in 54.5% of conditions, by 0.039%
  overall. So the mean is the right estimator and the minimum is not; the median of three
  carries about 1.35 times the variance of the mean for symmetric noise, which is 35% of the
  rounds thrown away. This repository has been reading medians.

The ten identical-code cells are: run 45594428603, which compares `efef397` with itself on both
machines and both programs; the x86-64 cells of runs 46075696888 and 46082844032, where the
command added is `[aarch64]` and x86-64 ran the same jar under the same flags; and the AArch64
cells of run 46039727701, where `BRANCH_FREE_MAX` was already true and the commit changes
nothing AArch64 compiles. Their spread is the noise floor, measured directly.

## The design constant

Writing `N` for the number of round-pairs behind one cell -- `k` rounds of each side in each of
`R` runs, so `N = k * R` -- the paired statistic has

    sd(D) = sigma_1 / sqrt(N),    sigma_1 = 2.02%

`sigma_1` is the worst cell, Neoverse-N2 full; the others are 1.73% to 1.86%. Because the
pooled CV is stable over time (1.32% over the earlier half of the runs, 1.12% over the later),
`sigma_1` can be taken from this history rather than re-estimated from three rounds, and the
interval uses the normal quantile instead of a t with two degrees of freedom. That alone is
worth 79% of the interval's width.

## How many rounds

Round-pairs needed per cell for the 99% interval to land inside `+-d`:

| margin `d` | power .80 | power .90 | power .974 |
|---|---|---|---|
| 0.50% | 243 | 291 | 377 |
| 0.75% | 108 | 129 | 168 |
| 1.00% | 61 | 73 | 95 |
| 1.25% | 39 | 47 | 61 |
| **1.50%** | 27 | **33** | **42** |
| 2.00% | 16 | 19 | 24 |
| 2.50% | 10 | 12 | 16 |
| 3.00% | 7 | 9 | 11 |

The last column is the one to use when the verdict is "flat everywhere", because that claim is
the conjunction of four cells. A conjunction of equivalence claims needs no correction to its
confidence -- each cell may be tested at the full 99%, and the combined claim still holds at
99%, by the intersection-union principle -- but its power is the product, so each cell must
reach .974 for the set to reach .90. The opposite claim, "a gain in at least one cell", is a
union and does need its alpha divided by four.

## What three rounds buys

| `N` | 99% interval | tightest flat claim | smallest detectable effect |
|---|---|---|---|
| **3 (current)** | **+-3.00%** | **\|d\| < 3.00%** | **3.98%** |
| 6 | +-2.12% | \|d\| < 2.12% | 2.82% |
| 12 | +-1.50% | \|d\| < 1.50% | 1.99% |
| 18 | +-1.23% | \|d\| < 1.23% | 1.63% |
| 24 | +-1.06% | \|d\| < 1.06% | 1.41% |
| 36 | +-0.87% | \|d\| < 0.87% | 1.15% |
| 48 | +-0.75% | \|d\| < 0.75% | 1.00% |

The present design is three rounds in one run, `N = 3`. Its 99% interval is 3 percentage points
wide either side, so the strongest flatness it can state is "smaller than 3%", and only in the
case where the measured difference comes out exactly zero. Its smallest detectable effect is
3.98%.

No experiment in this repository has ever measured an effect that large. Every verdict recorded
so far -- kept and reverted alike -- was read off an interval that contained both a substantial
gain and a substantial loss. That is the whole of the C series, and it is also `7f0c771` at
+0.26% and `cdd4447` at +0.51%, which were kept while measuring slower.

At `N = 3` the experiments do not have the power to support any verdict this project cares
about. Thirty-three round-pairs a cell, eleven times what is spent now, is what a 1.5% flat
verdict costs; forty-two if all four cells must carry it.

## Two ways to pay

**Buy the rounds.** 33 round-pairs is `k = 11` rounds in each of `R = 3` runs, or `k = 17` in
each of two. On the Neoverse-N2 full search, at 86 s a round, 33 round-pairs is 1.6 h of runner
time for that cell. Splitting across runs rather than lengthening one run is free -- the
run-level term cancels, so the only thing that matters is `k * R` -- and it keeps a single job
inside its time limit. It also lets the variance be checked against the between-run spread
instead of taken on faith.

**Or make the measurement quieter, which is cheaper.** `N` scales with the square of the noise,
and the noise is the JVM's: 1.2--1.6% against BLAST+'s 0.34% on the same host. Halving it to
0.6% would take `sigma_1` to about 0.85% and the 1.0% margin from 95 round-pairs to 17. Worth
testing, each against the round CV rather than against a reference commit: a fixed heap with
`-Xms` equal to `-Xmx` and `+AlwaysPreTouch`, so that heap growth and page faulting stop varying;
an AOT or AppCDS archive, so that the JIT does the same work each time; and pinning the worker
threads, since the Neoverse-N2 runner has four CPUs and the scheduler is free to move them.

Either way the first thing to schedule is a calibration run: one commit measured against itself,
as run 45594428603 did. It needs no reference, it costs one run, and it measures the noise floor
directly. Repeating it after any change to the harness is what keeps `sigma_1` honest.

## Rules that follow

- State the margin before the run, not after. "Flat" without a `d` is not a result.
- Fix `N` before the run. Adding rounds until the answer settles invalidates the interval.
- Report the mean of the rounds and its 99% interval. The median of three wastes a third of them.
- Report the interval, not just the point estimate. A verdict is the interval's position
  relative to `+-d`, and nothing else.
