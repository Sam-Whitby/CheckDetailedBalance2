#!/usr/bin/env wolframscript
(* ================================================================
   check.wls  --  Exhaustive detailed-balance checker
   ================================================================
   Usage:
     wolframscript -file check.wls <alg_file> [options]

   The .wl file must define:
     energy[state_]       bare energy function (no beta)
     Algorithm[state_]    MCMC move function using native random calls
     BitsToState[bits_]   {0,1,...} list -> state, or None
     numBeta              numeric inverse temperature

   Optional (for symbolic check with energy parameters):
     symParams = <|"eps" -> {...}, "couplings" -> {...}|>
       -- static symbolic parameters (always the same across components)

     DynamicSymParams[states_List] -> <|"couplings" -> {...}|>
       -- called once per component; return the J symbols needed for
          that component (e.g. auto-generated from particle types).
          Takes precedence over / merges with static symParams.

   Options (key=value, no spaces around =):
     MaxBitString=XXXX    Largest bit string to test (default 11111111)
                          All strings from length 1 through len(XXXX),
                          stopping at XXXX within the final length.
     Mode=Symbolic        Symbolic detailed-balance check only
     Mode=Numerical       Numerical MCMC check only
     Mode=Both            Both checks (default)
     NSteps=N             MCMC steps per seed state (default 50000)
     MaxBitDepth=N        BFS tree depth cap per state (default 20)
     Verbose=True/False   Print per-state BFS progress (default False)
   ================================================================ *)

scriptDir = DirectoryName[$InputFileName];
Get[scriptDir <> "dbc_core.wl"];

(* ---- Launch parallel kernels for FullSimplify parallelisation ---- *)
$nKernels = 0;
Quiet[Check[LaunchKernels[]; $nKernels = Length[Kernels[]], $nKernels = 0]];
If[$nKernels > 0,
  Print["Parallel: ", $nKernels, " kernels launched."],
  Print["Note: parallel kernel launch failed — running single-threaded."]];

(* ---- Parse command-line arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "check.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1,
  Print["Usage: wolframscript -file check.wls <alg_file> [options]"];
  Print[""];
  Print["Required definitions in the .wl file:"];
  Print["  energy[state_]       bare energy (no beta)"];
  Print["  Algorithm[state_]    MCMC move function"];
  Print["  BitsToState[bits_]   {0,1,...} -> state, or None if invalid"];
  Print["  numBeta              numeric inverse temperature"];
  Print[""];
  Print["Options (key=value, no spaces around =):"];
  Print["  MaxBitString=XXXX    Largest bit string (default 11111111)"];
  Print["  Mode=Symbolic|Numerical|Both  (default Both)"];
  Print["  NSteps=N             MCMC steps per system (default 50000)"];
  Print["  MaxBitDepth=N        BFS depth cap (default 20)"];
  Print["  Verbose=True/False   Per-state BFS output (default False)"];
  Exit[1]];

(* ---- Option defaults ---- *)
maxBitString = "11111111";
mode         = "Both";
nSteps       = 50000;
maxBitDepth  = 20;
verbose      = False;

Do[
  If[StringMatchQ[opt, "MaxBitString=" ~~ __], maxBitString = StringDrop[opt, 13]];
  If[StringMatchQ[opt, "Mode="         ~~ __], mode         =              StringDrop[opt,  5]];
  If[StringMatchQ[opt, "NSteps="       ~~ __], nSteps       = ToExpression[StringDrop[opt,  7]]];
  If[StringMatchQ[opt, "MaxBitDepth="  ~~ __], maxBitDepth  = ToExpression[StringDrop[opt, 12]]];
  If[StringMatchQ[opt, "Verbose="      ~~ __], verbose      = ToExpression[StringDrop[opt,  8]]],
  {opt, kvArgs}];

(* ---- Validate MaxBitString ---- *)
If[!StringMatchQ[maxBitString, ("0" | "1")..],
  Print["ERROR: MaxBitString must contain only 0s and 1s, got: ", maxBitString];
  Exit[1]];

(* ---- Resolve file path ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

(* ---- Load algorithm file ---- *)
ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState,
         fieldF, couplingJ,
         $fieldFConcrete, $couplingJConcrete, $abstractFunctions, $concreteParams,
         $fieldFreeParams, $couplingFreeParams];
Print["Loading: ", algFile];
Get[algFile];

(* ---- Validate required symbols ---- *)
If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[energy]] == 0,
  Print["ERROR: 'energy' not defined in ", posArgs[[1]]]; Exit[1]];
If[Length[DownValues[BitsToState]] == 0,
  Print["ERROR: 'BitsToState' not defined in ", posArgs[[1]]]; Exit[1]];
If[!ValueQ[numBeta] || !NumericQ[numBeta],
  Print["ERROR: 'numBeta' (numeric) not defined in ", posArgs[[1]]]; Exit[1]];

$alg         = Algorithm;
$energy      = energy;
$numBeta     = numBeta;
$hasDynParams = Length[DownValues[DynamicSymParams]] > 0;

(* ---- Verify couplingJ symmetry (if defined) ---- *)
(* Check type pairs {1..5} at d²∈{1,2,4,5,8,9,10} — covers all
   physically relevant cases.  Uses Simplify only: negligible overhead.
   An asymmetric coupling immediately implies a DB violation.
   In abstract-function mode ($abstractFunctions=True), couplingJ has no
   DownValues; we check $couplingJConcrete instead via a temporary Block. *)
$checkCouplingSymmetry[] := (
  $coupSymViolations = CheckCouplingSymmetry[{1,2,3,4,5}, {1,2,4,5,8,9,10}];
  If[$coupSymViolations =!= {},
    Print["ERROR: couplingJ is NOT symmetric — detailed balance cannot hold."];
    Print["  First violation: a=", $coupSymViolations[[1]]["a"],
          "  b=", $coupSymViolations[[1]]["b"],
          "  d2=", $coupSymViolations[[1]]["d2"],
          "  J[a,b]-J[b,a] = ", $coupSymViolations[[1]]["asymmetry"]];
    Exit[1]];
  Print["couplingJ symmetry: PASS"]);
Which[
  Length[DownValues[couplingJ]] > 0,
    $checkCouplingSymmetry[],
  Length[DownValues[$couplingJConcrete]] > 0,
    Block[{couplingJ},
      couplingJ[a_Integer, b_Integer, d2_Integer] := $couplingJConcrete[a, b, d2];
      $checkCouplingSymmetry[]]];

(* Static symbolic parameters (may be None if only DynamicSymParams is used) *)
$staticSysParams = If[ValueQ[symParams], symParams, None];
If[AssociationQ[$staticSysParams],
  $staticSymVars = Join[Lookup[$staticSysParams, "eps",       {}],
                        Lookup[$staticSysParams, "couplings", {}]],
  $staticSymVars = {}];

(* Fixed symbolic parameters: kept free during the symbolic check so that
   transition matrix entries contain exact symbolic expressions, but assigned
   their declared numeric values (not random) during numerical MCMC.
   Declared in the .wl file as:  symParams = <|"fixedParams" -> <|sym -> val, ...|>|> *)
$fixedParams  = If[AssociationQ[$staticSysParams],
                   Lookup[$staticSysParams, "fixedParams", <||>], <||>];
$fixedSymVars = Keys[$fixedParams];
$fixedNumVals = Values[$fixedParams];

(* ---- Generate candidate bit strings ---- *)
maxBitVec = ToExpression /@ Characters[maxBitString];
maxLen     = Length[maxBitVec];
maxVal     = FromDigits[maxBitVec, 2];

(* If the algorithm defines ValidStateIDs[maxId], use it to enumerate only
   the integer IDs that map to valid states (e.g. perfect-square arrays for
   the 2D lattice), skipping the bulk of bit strings that BitsToState would
   immediately reject.  Otherwise fall back to exhaustive enumeration. *)
allBitStrings =
  If[Length[DownValues[ValidStateIDs]] > 0,
    IntegerDigits[#, 2] & /@ ValidStateIDs[maxVal],
    Join[
      Flatten[Table[IntegerDigits[i, 2, len],
                    {len, 1, maxLen - 1}, {i, 0, 2^len - 1}], 1],
      Table[IntegerDigits[i, 2, maxLen], {i, 0, maxVal}]]];

(* ---- Column widths ---- *)
colBits  = Max[14, StringLength[maxBitString] + 2];
colState = 26; colN = 9; colSym = 10; colFail = 8; colKL = 12; colNum = 10;
totalCols = colBits + colState + colN +
  If[mode =!= "Numerical", colSym + colFail, 0] +
  If[mode =!= "Symbolic",  colKL  + colNum,  0];

(* ---- Print header ---- *)
Print[];
Print[StringRepeat["=", totalCols]];
Print["DETAILED BALANCE CHECKER"];
Print["File:          ", algFile];
Print["MaxBitString:  ", maxBitString,
      "  (", Length[allBitStrings], " bit strings generated)"];
Print["Mode:          ", mode];
If[mode =!= "Symbolic", Print["NSteps:        ", nSteps]];
Print["MaxBitDepth:   ", maxBitDepth];
Print[StringRepeat["=", totalCols]];
Print[];

Print[StringPadRight["Bits",    colBits],
      StringPadRight["State",   colState],
      StringPadRight["#States", colN],
      If[mode =!= "Numerical",
         StringPadRight["Symbolic", colSym] <> StringPadRight["#Fail", colFail], ""],
      If[mode =!= "Symbolic",
         StringPadRight["KL", colKL] <> StringPadRight["Numerical", colNum], ""]];
Print[StringRepeat["-", totalCols]];

(* ================================================================
   Progress bar helpers
   ================================================================ *)
$nBitStrings = Length[allBitStrings];
$checkStart  = AbsoluteTime[];
$barWidth    = 28;

formatDuration[secs_] :=
  Which[secs < 60,   ToString[Round[secs]] <> "s",
        secs < 3600, ToString[Floor[secs/60]] <> "m " <>
                     ToString[Round[Mod[secs, 60]]] <> "s",
        True,        ToString[Floor[secs/3600]] <> "h " <>
                     ToString[Floor[Mod[secs/60, 60]]] <> "m"];

writeProgress[idx_] := Module[{elapsed, pct, filled, arrow, bar, eta, etaStr, elStr},
  elapsed = AbsoluteTime[] - $checkStart;
  pct     = If[$nBitStrings > 0, idx / $nBitStrings, 0.];
  filled  = Round[pct * $barWidth];
  arrow   = If[filled < $barWidth, ">", ""];
  bar     = "[" <> StringRepeat["=", filled] <> arrow <>
            StringRepeat[" ", Max[0, $barWidth - filled - StringLength[arrow]]] <> "]";
  eta     = If[idx > 0 && elapsed > 0,
               elapsed * ($nBitStrings - idx) / idx,
               Infinity];
  etaStr  = If[!NumberQ[eta] || eta === Infinity, "?", formatDuration[eta]];
  elStr   = formatDuration[elapsed];
  WriteString[$Output,
    "\r" <> bar <> " " <> ToString[idx] <> "/" <> ToString[$nBitStrings] <>
    "  " <> ToString[Round[pct * 100]] <> "%  elapsed: " <> elStr <>
    "  ETA: " <> etaStr <> "   "]  (* trailing spaces clear previous longer text *)
];

(* Before printing a permanent row: clear the progress line first *)
clearProgress[] := WriteString[$Output, "\r" <> StringRepeat[" ", 78] <> "\r"];

(* ================================================================
   Main loop
   ================================================================ *)
coveredAssoc = <||>;   (* Association for O(1) covered-state lookup *)
testedCount  = 0;
symPassCount  = 0; symFailCount  = 0;
numPassCount  = 0; numFailCount  = 0;
$loopIdx     = 0;

(* ---- BFS / FS pipeline state ----
   After each component's BFS, FS tasks are submitted to subkernels
   asynchronously.  The next component's BFS then runs on the main
   kernel while FS evaluates in parallel.  flushPending[] collects
   the FS result (WaitAll) and prints it, called right after the
   following component's BFS so FS has had BFS-time to evaluate. *)
$fsPending = None;   (* record from SubmitDetailedBalanceFSTasks, or None *)
$fsMeta    = None;   (* print metadata for the pending component *)

flushPending[] := Module[{violations, symResult, failStr},
  If[$fsPending =!= None,
    violations = CollectDetailedBalanceFSTasks[$fsPending];
    symResult  = If[violations === {}, "PASS", "FAIL"];
    failStr    = If[violations === {}, "0", ToString[Length[violations]]];
    If[symResult === "PASS", symPassCount++, symFailCount++];
    clearProgress[];
    Print[StringPadRight[$fsMeta["bitsStr"],           colBits],
          StringPadRight[$fsMeta["stateStr"],          colState],
          StringPadRight[ToString[$fsMeta["nStates"]], colN],
          If[mode =!= "Numerical",
             StringPadRight[symResult, colSym] <> StringPadRight[failStr, colFail], ""],
          If[mode =!= "Symbolic",
             StringPadRight[$fsMeta["klStr"],     colKL] <>
             StringPadRight[$fsMeta["numResult"], colNum], ""]];
    $fsPending = None;
    $fsMeta    = None]];

Do[
  $loopIdx++;
  state = Check[BitsToState[bitString], $Failed];
  If[state === $Failed || state === None,
    writeProgress[$loopIdx]; Continue[]];
  If[KeyExistsQ[coveredAssoc, state],
    writeProgress[$loopIdx]; Continue[]];

  bitsStr  = StringJoin[ToString /@ bitString];
  stateStr = If[Length[DownValues[DisplayState]] > 0,
               DisplayState[state], ToString[state]];

  (* ---- BFS: discover connected component ---- *)
  treeData = BuildTreeAT[state, $alg,
               "MaxBitDepth" -> maxBitDepth,
               "TimeLimit"   -> 120.,
               "Verbose"     -> verbose];

  (* Collect and print the PREVIOUS component's FS result now that BFS
     is complete — FS has had BFS-time to evaluate before WaitAll blocks. *)
  flushPending[];

  If[MatchQ[treeData, $dbc$cantHandle[_]],
    clearProgress[];
    Print[StringPadRight[bitsStr,  colBits], StringPadRight[stateStr, colState],
          StringPadRight["?", colN], "ERROR: ", treeData[[1]]];
    writeProgress[$loopIdx];
    testedCount++; Continue[]];

  allStates = Sort @ Keys[treeData];
  matrix    = TreeATToMatrix[treeData];
  nStates   = Length[allStates];
  coveredAssoc = Join[coveredAssoc, Association[# -> True & /@ allStates]];
  testedCount++;

  (* ---- Per-component symbolic parameters ---- *)
  (* If DynamicSymParams is defined, call it with this component's states
     to get exactly the symbolic variables needed.
     "couplings" / "eps" — symbolic atoms for FullSimplify assumptions.
       For standard files these are scalar symbols like Jpair12.
       For abstract-function files these are function-application atoms
       like fieldF[0,0,2] and couplingJ[1,2,1]; FullSimplify treats each
       unevaluated application as an independent real atom.
     "numericParams" — scalar symbols for the numerical MCMC Block.
       For standard files same as "couplings".
       For abstract-function files: the $jPairSym symbols plus the
       concrete-implementation parameters (fieldAmp, lambdaJ, etc.). *)
  If[$hasDynParams,
    dynP         = DynamicSymParams[allStates];
    dynVars      = Join[Lookup[dynP, "eps", {}], Lookup[dynP, "couplings", {}]];
    dynNumVars   = Lookup[dynP, "numericParams", dynVars];
    compSymVars    = Union[$staticSymVars, dynVars];
    compNumericVars = Union[$staticSymVars, dynNumVars],
    compSymVars    = $staticSymVars;
    compNumericVars = $staticSymVars];
  (* Fixed params are real for FullSimplify purposes (e.g. jumpSigma ∈ Reals) *)
  compSymAssumptions = Map[Element[#, Reals] &, Union[compSymVars, $fixedSymVars]];

  (* ---- Submit FS asynchronously; MCMC runs while subkernels evaluate ---- *)
  $fsPending = If[mode =!= "Numerical",
    SubmitDetailedBalanceFSTasks[matrix, allStates, $energy, compSymAssumptions],
    None];

  (* ---- Numerical MCMC check (overlaps with async FS on subkernels) ---- *)
  klStr = "-"; numResult = "-";
  If[mode =!= "Symbolic" && nStates > 1,
    (* compNumericVars contains the scalar symbols that need numeric values:
       for standard files these are the coupling symbols (Jpair12 etc.);
       for abstract-function files these are $jPairSym symbols plus the
       concrete-implementation parameters (fieldAmp, lambdaJ, etc.). *)
    allBlockVars = Union[compNumericVars, $fixedSymVars];
    (* Helper that runs the MCMC inside the current Block scope *)
    $runNumericalMCMC[] := (
      mcmcCounts = RunNumericalMCMCAT[allStates, $alg, $numBeta, "NSteps" -> nSteps];
      boltzmann  = BoltzmannWeightsAT[allStates, $energy, $numBeta]);
    (* For abstract-function files: activate concrete fieldF / couplingJ
       implementations inside an outer Block before running MCMC. *)
    $activateAndRun[] :=
      If[ValueQ[$abstractFunctions] && TrueQ[$abstractFunctions],
        Block[{fieldF, couplingJ},
          If[Length[DownValues[$fieldFConcrete]] > 0,
            fieldF[x_Integer, y_Integer, L_Integer] := $fieldFConcrete[x, y, L]];
          If[Length[DownValues[$couplingJConcrete]] > 0,
            couplingJ[a_Integer, b_Integer, d2_Integer] := $couplingJConcrete[a, b, d2]];
          $runNumericalMCMC[]],
        $runNumericalMCMC[]];
    If[allBlockVars =!= {},
      SeedRandom[42];
      compNumVals = If[compNumericVars =!= {}, RandomReal[{-1, 1}, Length[compNumericVars]], {}];
      Block[Evaluate[allBlockVars],
        If[compNumericVars =!= {}, MapThread[Set, {compNumericVars, compNumVals}]];
        If[$fixedSymVars =!= {}, MapThread[Set, {$fixedSymVars, $fixedNumVals}]];
        $activateAndRun[]],
      (* Fully numeric: no symbolic vars of any kind *)
      $activateAndRun[]];
    simFreq = N[# / Total[mcmcCounts]] & /@ mcmcCounts;
    kl = Total @ Table[
      With[{p = simFreq[s], q = N @ boltzmann[s]},
        If[p > 0 && q > 0, p * Log[p / q], 0.]], {s, allStates}];
    klStr     = Which[
                  Abs[kl] < 5.*^-4, "~0",
                  True,              ToString[NumberForm[kl, {4, 3}]]];
    numResult = If[kl < 0.02, "PASS", "FAIL"];
    If[numResult === "PASS", numPassCount++, numFailCount++]];

  (* ---- Store metadata for deferred printing by flushPending[] ---- *)
  $fsMeta = <|"bitsStr"  -> bitsStr,  "stateStr" -> stateStr,
               "nStates"  -> nStates,
               "klStr"    -> klStr,   "numResult" -> numResult|>;

  (* Numerical-only mode: no FS submitted, print result row immediately *)
  If[$fsPending === None,
    clearProgress[];
    Print[StringPadRight[bitsStr,           colBits],
          StringPadRight[stateStr,          colState],
          StringPadRight[ToString[nStates], colN],
          If[mode =!= "Symbolic",
             StringPadRight[klStr, colKL] <> StringPadRight[numResult, colNum], ""]];
    $fsMeta = None];

  writeProgress[$loopIdx],

  {bitString, allBitStrings}];

(* Collect and print the last component's deferred FS result *)
flushPending[];

(* Final newline to clear progress bar before the summary *)
WriteString[$Output, "\n"];

(* ================================================================
   Summary
   ================================================================ *)
Print[];
Print[StringRepeat["-", totalCols]];
Print["Systems tested (distinct connected components): ", testedCount];
If[mode =!= "Numerical",
  Print["Symbolic  -- PASS: ", symPassCount, "  FAIL: ", symFailCount]];
If[mode =!= "Symbolic" && (numPassCount + numFailCount) > 0,
  Print["Numerical -- PASS: ", numPassCount, "  FAIL: ", numFailCount]];
If[symFailCount == 0 && numFailCount == 0 && testedCount > 0,
  Print["All tested systems satisfy detailed balance."],
  If[testedCount == 0,
    Print["No valid states found for the given MaxBitString."],
    Print["*** DETAILED BALANCE VIOLATION DETECTED ***"]]];
Print[StringRepeat["=", totalCols]];
