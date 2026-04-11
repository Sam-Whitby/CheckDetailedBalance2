(* ================================================================
   1D Kawasaki with wrong-sign acceptance -- FAILS detailed balance
   ================================================================

   Identical to kawasaki_1d.wl except for one line in Algorithm:

     CORRECT:  dE = energy[newState] - energy[state]
     BUGGY:    dE = energy[state] - energy[newState]

   With the wrong sign, MetropolisProb receives -ΔE instead of ΔE:
   uphill moves (to higher energy) are always accepted (prob 1), while
   downhill moves are penalised with exp(-β*|ΔE|). This inverts the
   Boltzmann distribution: the algorithm preferentially visits
   high-energy states rather than low-energy ones.

   For any two states with different energies E(i) ≠ E(j):
     T(i→j)·π(i) = P · exp(-β·E(i))   [always accepts uphill i→j]
     T(j→i)·π(j) = P · exp(-β·ΔE) · exp(-β·E(j)) = P · exp(-2β·E(j))
   These are equal only if E(i) = E(j), so detailed balance is violated
   for all state pairs with different energies.

   Note: a biased bond-selection probability (e.g. proposing bond 0
   more often) does NOT break detailed balance for Kawasaki, because
   the forward and reverse transitions of the same bond always use the
   same proposal probability, which cancels in the DB condition.
   ================================================================ *)


(* ---- Encoding (identical to kawasaki_1d.wl) ------------------------------ *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n - i]; elems = Delete[elems, idx + 1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i - 1]]; r = Mod[r, Factorial[i - 1]];
       AppendTo[perm, elems[[idx + 1]]]; elems = Delete[elems, idx + 1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}]; arr]

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[couplingJ[state[[i]], state[[Mod[i, L] + 1]], 1], {i, L}]]]


(* ---- Buggy algorithm: dE computed with wrong sign ----------------------- *)

Algorithm[state_List] :=
  Module[{L, b, s1, s2, newState, dE},
    L  = Length[state];
    b  = RandomInteger[{0, L - 1}];
    s1 = b + 1; s2 = Mod[b + 1, L] + 1;
    newState = ReplacePart[state, {s1 -> state[[s2]], s2 -> state[[s1]]}];
    dE = energy[state] - energy[newState];  (* BUG: sign reversed -- prefers high energy *)
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1
