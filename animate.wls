#!/usr/bin/env wolframscript
(* ================================================================
   animate.wls  --  Run an algorithm and animate it via Python/matplotlib
   ================================================================
   Usage:
     wolframscript -file animate.wls <alg_file> \
       Sites=<n> N=<n> [Steps=200] [Beta=1] [FPS=10] \
       [Simple=1] [RecordEvery=1] [J12=0.5 ...]

   Runs the algorithm for the specified number of steps, collects
   the state and energy at each recorded step, writes the data to
   a temporary JSON file, then launches animate_plot.py to display
   an animated matplotlib figure with:
     - Left panel:  imshow of the lattice (colour = particle type)
     - Right panel: energy vs step (growing line)

   For 2D algorithms (those defining DisplayState), the lattice is
   shown as an L x L grid. For 1D, it is shown as a 1 x Sites row.

   Arguments:
     Sites=<n>          Total lattice sites (e.g. 6 for a 1D ring,
                        9 for a 3x3 grid).
     N=<n>              Number of labeled particles (types 1..N).
     Steps=<n>          MCMC steps to run (default 200).
     Beta=<f>           Inverse temperature (overrides numBeta in file).
     FPS=<f>            Animation frame rate in frames per second (default 10).
     Simple=1           Use fast 2-colour mode: holes vs particles.
                        Renders much faster, enabling high frame rates.
     RecordEvery=<n>    Record state every n steps (default 1).
     J<a><b>=<f>        Set coupling constant J_ab (random if omitted).

   Examples:
     wolframscript -file animate.wls examples3/kawasaki_1d.wl \
       Sites=6 N=3 Steps=500 FPS=30 Simple=1

     wolframscript -file animate.wls examples3/kawasaki_2d.wl \
       Sites=9 N=4 Steps=300 FPS=10 J12=1.5

     wolframscript -file animate.wls examples3/cluster_2d.wl \
       Sites=9 N=3 Steps=400 FPS=15 RecordEvery=2
   ================================================================ *)

(* Real MetropolisProb for actual execution (not BFS interception) *)
MetropolisProb[dE_] := If[dE <= 0, 1., Exp[-$animBeta * dE]]

scriptDir = DirectoryName[$InputFileName];

(* ---- Parse arguments ---- *)
args = $ScriptCommandLine;
If[Length[args] >= 1 && StringMatchQ[args[[1]], ___ ~~ "animate.wls"],
  args = Rest[args]];

posArgs = Select[args, !StringMatchQ[#, __ ~~ "=" ~~ __] &];
kvArgs  = Select[args,  StringMatchQ[#, __ ~~ "=" ~~ __] &];

If[Length[posArgs] < 1 ||
   !MemberQ[kvArgs, _?(StringMatchQ[#, "Sites=" ~~ __] &)] ||
   !MemberQ[kvArgs, _?(StringMatchQ[#, "N=" ~~ __] &)],
  Print["Usage: wolframscript -file animate.wls <alg_file> Sites=<n> N=<n> [options]"];
  Print["  Sites=<n>        total lattice sites"];
  Print["  N=<n>            number of labeled particles (types 1..N)"];
  Print["  Steps=<n>        MCMC steps (default 200)"];
  Print["  Beta=<f>         inverse temperature (default from .wl file)"];
  Print["  FPS=<f>          animation frame rate in frames per second (default 10)"];
  Print["  Simple=1         fast 2-colour mode: holes vs particles"];
  Print["  RecordEvery=<n>  record every nth step (default 1)"];
  Print["  J<a><b>=<f>      coupling constant (random if omitted)"];
  Exit[1]];

(* ---- Defaults ---- *)
nSites      = 4;
nPart       = 2;
nSteps      = 200;
animBeta    = None;
fps         = 10.;
simpleMode  = False;
recordEvery = 1;

Do[
  If[StringMatchQ[opt, "Sites="       ~~ __], nSites      = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "N="           ~~ __], nPart       = ToExpression[StringDrop[opt,  2]]];
  If[StringMatchQ[opt, "Steps="       ~~ __], nSteps      = ToExpression[StringDrop[opt,  6]]];
  If[StringMatchQ[opt, "Beta="        ~~ __], animBeta    = ToExpression[StringDrop[opt,  5]]];
  If[StringMatchQ[opt, "FPS="         ~~ __], fps         = N[ToExpression[StringDrop[opt,  4]]]];
  (* Legacy Delay= support: convert seconds -> fps *)
  If[StringMatchQ[opt, "Delay="       ~~ __], fps         = N[1 / ToExpression[StringDrop[opt, 6]]]];
  If[StringMatchQ[opt, "Simple="      ~~ __], simpleMode  = (StringDrop[opt, 7] =!= "0")];
  If[StringMatchQ[opt, "RecordEvery=" ~~ __], recordEvery = ToExpression[StringDrop[opt, 12]]],
  {opt, kvArgs}];

(* ---- Load algorithm file ---- *)
algFile = posArgs[[1]];
If[!FileExistsQ[algFile], algFile = scriptDir <> algFile];
If[!FileExistsQ[algFile],
  Print["ERROR: File not found: ", posArgs[[1]]]; Exit[1]];

ClearAll[Algorithm, energy, BitsToState, numBeta, symParams, DynamicSymParams, DisplayState];
Get[algFile];

If[Length[DownValues[Algorithm]] == 0,
  Print["ERROR: 'Algorithm' not defined in ", algFile]; Exit[1]];

$animBeta = If[animBeta =!= None, animBeta, numBeta];

(* ---- Build initial state ---- *)
If[nPart > nSites,
  Print["ERROR: N=", nPart, " exceeds Sites=", nSites]; Exit[1]];
SeedRandom[1234];
initialState = RandomSample[Join[Range[nPart], ConstantArray[0, nSites - nPart]], nSites];

(* ---- Determine display grid shape ---- *)
(* 2D if DisplayState is defined (convention for 2D algorithm files) *)
is2D   = Length[DownValues[DisplayState]] > 0;
gridL  = If[is2D, Round[Sqrt[nSites]], nSites];
nRows  = If[is2D, gridL, 1];
nCols  = If[is2D, gridL, nSites];

(* ---- Assign coupling constants ---- *)
(* Parse J<a><b>=<f> args *)
jSupplied = Association @ Map[
  Function[opt, With[{p = StringSplit[opt, "=", 2]},
    ToExpression[p[[1]]] -> ToExpression[p[[2]]]]],
  Select[kvArgs, StringMatchQ[#, "J" ~~ DigitCharacter ~~ DigitCharacter ~~ "=" ~~ __] ||
                 StringMatchQ[#, "J" ~~ DigitCharacter ~~ "=" ~~ __] &]];

types = Sort[DeleteCases[initialState, 0]];
(* Random J values drawn from N(0,1) so bonds can be attractive or repulsive *)
SeedRandom[42];
jAssigned = <||>;
Do[If[a < b,
  sym = ToExpression["J" <> ToString[a] <> ToString[b]];
  val = Lookup[jSupplied, sym, RandomVariate[NormalDistribution[0, 1]]];
  Set[Evaluate[sym], val];
  AssociateTo[jAssigned, ToString[sym] -> val]],
  {a, types}, {b, types}];

(* ---- Run algorithm and collect data ---- *)
Print["Running ", nSteps, " steps of ", algFile, " ..."];
state = initialState;
stepsRec    = {0};
statesRec   = {Normal[state]};      (* ensure plain list, not packed *)
energiesRec = {N[energy[state]]};

Do[
  state = Algorithm[state];
  If[Mod[i, recordEvery] == 0,
    AppendTo[stepsRec,    i];
    AppendTo[statesRec,   Normal[state]];
    AppendTo[energiesRec, N[energy[state]]]],
  {i, 1, nSteps}];

Print["Recorded ", Length[stepsRec], " frames. Launching plot..."];

(* ---- Write JSON data file ---- *)
dataFile = FileNameJoin[{$TemporaryDirectory, "dbc_animate_data.json"}];
Export[dataFile,
  <|"algo_file"  -> algFile,
    "grid_rows"  -> nRows,
    "grid_cols"  -> nCols,
    "n_types"    -> If[types === {}, 0, Max[types]],
    "fps"        -> fps,
    "simple"     -> simpleMode,
    "steps"      -> stepsRec,
    "states"     -> statesRec,
    "energies"   -> energiesRec,
    "params"     -> Join[<|"beta" -> $animBeta|>, jAssigned]|>];

(* ---- Launch Python animation ---- *)
pyScript = FileNameJoin[{scriptDir, "animate_plot.py"}];
If[!FileExistsQ[pyScript],
  Print["ERROR: animate_plot.py not found at ", pyScript]; Exit[1]];

result = RunProcess[{"python3", pyScript, dataFile}];
If[result["ExitCode"] != 0,
  Print["Python error:\n", result["StandardError"]];
  Print["Standard output:\n", result["StandardOutput"]]];
