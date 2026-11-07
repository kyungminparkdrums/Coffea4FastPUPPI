import os, glob, yaml, argparse, importlib
import multiprocessing as mp
import faulthandler
from coffea import processor
from coffea.nanoevents import NanoAODSchema
import coffea.util as util

from processor import P2L1TAnalyzer
from histo_config.expand_histo import expand_histo_yaml


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rootDir", required=True, help="Directory OR single .root file")
    parser.add_argument("--filePattern", default="*.root", help="Used if rootDir is a directory")
    parser.add_argument("--treename", default="Events")

    parser.add_argument("--outName", required=True)
    parser.add_argument("--histoConfig", required=True)

    parser.add_argument("--cutConfig", required=True,
                        help="Python module path, e.g. cut_config.cut_config")

    parser.add_argument("--workers", type=int, default=4)

    parser.add_argument("--executor", choices=["iterative", "futures"], default="iterative",
                        help="Serial by default for ML inference diagnosis")
    parser.add_argument("--chunksize", type=int, default=10000)
    parser.add_argument("--maxchunks", type=int, default=None)
    parser.add_argument("--maxfiles", type=int, default=500)
    parser.add_argument("--merging", action="store_true",
                        help="Enable parallel merging with the futures executor")

    args = parser.parse_args()
    for name in ("workers", "chunksize", "maxchunks", "maxfiles"):
        value = getattr(args, name)
        if value is not None and value < 1:
            parser.error(f"--{name} must be positive")
    faulthandler.enable(all_threads=True)
    if args.executor == "futures":
        mp.set_start_method("spawn", force=True)


    # hist yaml
    with open(args.histoConfig) as f:
        raw = yaml.safe_load(f)
    hist_config = expand_histo_yaml(raw)

    # load cut config module
    cutflow_module = importlib.import_module(args.cutConfig)

    # files
    if args.rootDir.endswith(".root") and os.path.isfile(args.rootDir):
        file_list = [args.rootDir]
    else:
        file_list = sorted(glob.glob(os.path.join(args.rootDir, args.filePattern)))
        file_list = [f for f in file_list if "perfNano_16420723_143.root" not in f]
        file_list = file_list[:args.maxfiles]

    if not file_list:
        raise RuntimeError(f"No ROOT files found: rootDir={args.rootDir}, pattern={args.filePattern}")

    fileset = {"files": {"files": file_list}}

    print(f"Input: {args.rootDir} ({len(file_list)} files), treename={args.treename}")
    print(f"Stages: {[name for name, _ in cutflow_module.CUTFLOW]}")

    proc = P2L1TAnalyzer(hist_config=hist_config, cut_config_module=cutflow_module)

    if args.executor == "iterative":
        exe = processor.IterativeExecutor()
    else:
        exe = processor.FuturesExecutor(workers=args.workers, merging=args.merging)
    print(f"Executor={args.executor}, workers={args.workers if args.executor == 'futures' else 1}, "
          f"chunksize={args.chunksize}, maxchunks={args.maxchunks}, "
          f"ML device={getattr(cutflow_module, 'ML_DEVICE', 'not specified')}, "
          f"ML batch={getattr(cutflow_module, 'ML_INFER_BATCH_GRAPHS', 'not specified')}")
    runner = processor.Runner(executor=exe, schema=NanoAODSchema,
                              chunksize=args.chunksize, maxchunks=args.maxchunks)

    result = runner(fileset=fileset, treename=args.treename, processor_instance=proc)
    
    # ---- build plot_cfg sidecar  ----
    plot_cfg = {}
    hists_by_stage = result.get("hists", {}) or {}

    for stage, hdict in hists_by_stage.items():
        plot_cfg[stage] = {}
        for hname in hdict.keys():
            cfg = hist_config.get(hname, {}) or {}
            plot_cfg[stage][hname] = {"logy": bool(cfg.get("logy", False))}

    result["plot_cfg"] = plot_cfg

    os.makedirs("coffea", exist_ok=True)
    outpath = f"coffea/{args.outName}.coffea"
    util.save(result, outpath)
    print("Saved:", outpath)


if __name__ == "__main__":
    main()

