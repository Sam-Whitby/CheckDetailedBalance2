(* ================================================================
   1D jump dynamics on a periodic ring -- PASSES detailed balance
   ================================================================

   State:   flat array of length L.
            State[[i]] = 0 (empty) or k \[Element] {1,...,N} (labeled particle).
            L and N are inferred from the bit string by BitsToState.

   Encoding: bit string \[Rule] integer ID \[Rule] unique labeled array.
             ID encodes (L, N, particle positions, label permutation).
             See $decode below.

   Move:    Pick a particle uniformly; propose a displacement d drawn
            uniformly from {0,...,L-1}; accept the displacement with
            probability proportional to the normal bin weight
              w(d) = [Erf((sd+1/2)/(sqrt(2) sigma)) - Erf((sd-1/2)/(sqrt(2) sigma))] / 2
            where sd is the signed displacement on the ring
            (sd = d if 2d <= L, else sd = d - L).
            This is the probability that Normal(0, sigma^2) falls in
            the bin of width 1 centred on sd, so nearby sites are
            preferred over distant ones following the normal profile.
            sigma = jumpSigma (set at the bottom of this file).
            The proposal is symmetric: w(d) = w(L-d) by the symmetry
            of Erf, so the Metropolis acceptance ratio reduces to
            exp(-beta*dE) exactly.
            If the target site is occupied the move is hard-rejected.

   Energy:  Nearest-neighbour pairwise coupling J_ab between particle
            types a < b. Holes (type 0) do not contribute to energy.

   Implementation note for the detailed-balance checker:
            Displacement sampling is done as rejection sampling:
            RandomInteger[{0,L-1}] selects d uniformly, then
            RandomReal[] >= w(d)/w(0) rejects with probability
            1 - w(d)/w(0).  Both primitives are intercepted correctly
            by RunWithBitsAT (RandomInteger via bit-read rejection
            sampling; RandomReal via the acceptTest UpValues).
            RandomChoice[weights->elements] is NOT used because the
            checker treats a Rule as a 2-element list and returns
            the weights or elements list rather than a chosen element.
   ================================================================ *)


(* ---- Bijective integer encoding ----------------------------------------- *)
(* Maps every labeled-particle array of any length to a unique integer.
   ID = countPrefix[L] + countNPrefix[L,N] + rankCombo * N! + rankPerm
   where countPrefix accumulates array counts over all smaller lengths. *)

$cL[L_]       := $cL[L]    = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]    := $cLPre[L] = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_]:= $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

(* Combinatorial number system: rank/unrank sorted 0-indexed N-combinations *)
$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

(* Factorial number system: rank/unrank permutations of {1,...,N} *)
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


(* ---- Coupling function --------------------------------------------------- *)

couplingJ[a_Integer, b_Integer, d2_Integer] :=
  If[a == 0 || b == 0 || d2 != 1, 0, $jPairSym[a, b]]

DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]]},
    <|"couplings" ->
      Flatten @ Table[
        If[a < b, {$jPairSym[a, b]}, Nothing],
        {a, types}, {b, types}]|>]


(* ---- Energy -------------------------------------------------------------- *)

energy[state_List] :=
  With[{L = Length[state]},
    Total[Table[couplingJ[state[[i]], state[[Mod[i, L] + 1]], 1], {i, L}]]]


(* ---- Normal bin weight --------------------------------------------------- *)

(* Probability that Normal(0, sigma^2) falls in the unit bin centred on the
   signed displacement sd corresponding to index k on a ring of size L.
   Signed displacement: sd = k if 2k <= L, else sd = k - L.
   This is the un-normalised proposal weight for a jump of k sites.
   w(k) = w(L-k) by the antisymmetry of Erf, so the proposal is symmetric. *)
$normalBinW[k_Integer, sigma_, L_Integer] :=
  With[{sd = If[2 k <= L, k, k - L]},
    (Erf[(sd + 1/2) / (Sqrt[2] sigma)] - Erf[(sd - 1/2) / (Sqrt[2] sigma)]) / 2]


(* ---- Algorithm ----------------------------------------------------------- *)

(* Rejection-sampling implementation of a discretised Normal jump.
   RandomInteger[{0,L-1}] proposes d uniformly; the displacement is
   then thinned by RandomReal[] >= w(d)/w(0).  Both calls are handled
   correctly by the checker's RunWithBitsAT interception mechanism.
   The occupancy check comes first so that d=0 (always occupied) never
   consumes a RandomReal bit, keeping the BFS tree compact. *)
Algorithm[state_List] :=
  Module[{L, particlePositions, nParticles, pickedIdx, pickedPos,
          d, newPos, newState, dE},
    L = Length[state];
    particlePositions = Flatten[Position[state, _?(# != 0 &)]];
    nParticles = Length[particlePositions];
    If[nParticles == 0, Return[state]];

    (* Choose a particle uniformly at random *)
    pickedIdx = RandomInteger[{1, nParticles}];
    pickedPos = particlePositions[[pickedIdx]];

    (* Choose a trial displacement uniformly from {0, 1, ..., L-1} *)
    d = RandomInteger[{0, L - 1}];

    (* Target site under periodic boundary conditions (1-indexed) *)
    newPos = Mod[pickedPos - 1 + d, L] + 1;

    (* Hard rejection if target is occupied (catches d=0 without extra RNG) *)
    If[state[[newPos]] != 0, Return[state]];

    (* Thin the uniform proposal to match the normal distribution.
       w(0) = $normalBinW[0,...] is the maximum weight (normal peaks at 0),
       so the acceptance ratio w(d)/w(0) is always in (0, 1].
       Symmetry w(d) = w(L-d) means the forward and reverse acceptance
       probabilities are equal, so the proposal ratio = 1. *)
    If[RandomReal[] >= $normalBinW[d, jumpSigma, L] / $normalBinW[0, jumpSigma, L],
       Return[state]];

    (* Build candidate state: particle moves from pickedPos to newPos *)
    newState = ReplacePart[state, {pickedPos -> 0, newPos -> state[[pickedPos]]}];

    (* Metropolis acceptance *)
    dE = energy[newState] - energy[state];
    If[RandomReal[] < MetropolisProb[dE], newState, state]]


(* ---- Checker interface --------------------------------------------------- *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2]},
    If[id == 0, None, $decode[id]]]

numBeta = 1

(* Standard deviation of the jump distribution in lattice units.
   Kept as a free symbol during the symbolic check so that the transition
   matrix entries display as exact Erf expressions.  Assigned the value 2
   (in lattice units) during numerical MCMC and animations via the
   "fixedParams" mechanism in symParams.
   The sigma-dependent factors cancel in the detailed-balance ratio
   because w(d) = w(L-d) holds as a structural algebraic identity. *)
symParams = <|"fixedParams" -> <|jumpSigma -> 2|>|>
