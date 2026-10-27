"""One public command for inference, explanation and independent evaluation."""

import argparse
import json
import logging
import os
from pathlib import Path

from mendelprime.domain import MendelPrimeError

from .config import apply_overrides, load_config
from .io import read_patients
from .runner import Runner, evaluate, result_paths

LOG = logging.getLogger(__name__)


def parser():
    root = argparse.ArgumentParser(prog="mendelprime")
    commands = root.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("run", "rank canonical case inputs"),
        ("explain", "generate explanations without reranking"),
    ):
        command = commands.add_parser(name, help=help_text)
        for option in ("input", "config", "output-root"):
            command.add_argument(f"--{option}", type=Path, required=True)
        command.add_argument("--cache-root", type=Path)
        command.add_argument("--audit-root", type=Path, help="Optional private model-call audit directory")
        command.add_argument("--cohort")
        command.add_argument("--lanes", type=int)
        command.add_argument("--max-llm-concurrency", type=int)
        if name == "run":
            command.add_argument("--explanation", action=argparse.BooleanOptionalAction, default=None)
        else:
            command.add_argument(
                "--case-inputs", type=Path, help="Original canonical inputs for fresh family-wide QC"
            )
            command.add_argument(
                "--recompute-assessment",
                action="store_true",
                help="rebuild post-ranking ACMG evidence and report; never rerun rankings",
            )
            command.add_argument(
                "--rerender-report",
                action="store_true",
                help="regenerate a saved report without rerunning ranking or a saved assessment",
            )
    evaluation = commands.add_parser("evaluate", help="compare public rankings with external Gold labels")
    evaluation.add_argument("--input", type=Path, required=True)
    evaluation.add_argument("--truth", type=Path, required=True)
    evaluation.add_argument("--cohort")
    diagnosis = commands.add_parser("doctor", help="check configured resources without patient or model calls")
    diagnosis.add_argument("--config", type=Path, required=True)
    diagnosis.add_argument("--profile", choices=("hpo", "vcf"), required=True)
    diagnosis.add_argument("--explanation", action=argparse.BooleanOptionalAction, default=None)
    diagnosis.add_argument("--check-tools", action="store_true", help="also launch tool version/import checks")
    diagnosis.add_argument("--work-root", type=Path, help="Existing directory visible inside a configured VEP container")
    diagnosis.add_argument("--portable-root", type=Path, help="Reject resource paths or symlinks escaping this bundle")
    return root


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if os.environ.get("MENDELPRIME_VERIFY_TLS") == "0":
        LOG.warning("HTTPS certificate verification disabled by MENDELPRIME_VERIFY_TLS=0")
    try:
        args = parser().parse_args(argv)
    except SystemExit as error:
        return int(error.code or 0)
    except KeyboardInterrupt:
        return 130
    try:
        if args.command == "doctor":
            from .deployment import doctor

            config = apply_overrides(load_config(args.config), explanation=args.explanation)
            summary = doctor(config, args.profile, check_tools=args.check_tools,
                             work_root=args.work_root, portable_root=args.portable_root)
            code = summary["exit_code"]
        elif args.command == "evaluate":
            summary = evaluate(result_paths(args.input), args.truth, args.cohort)
            code = 0
        else:
            config = apply_overrides(
                load_config(args.config),
                lanes=args.lanes,
                max_llm_concurrency=args.max_llm_concurrency,
                explanation=getattr(args, "explanation", None),
            )
            cache = (args.cache_root or args.output_root.parent / ".cache").absolute()
            runner = Runner(config, cache, audit_root=args.audit_root)
            if args.command == "run":
                patients = read_patients(args.input, args.cohort)
                summary = runner.run(patients, args.output_root)
            else:
                summary = runner.explain(
                    result_paths(args.input),
                    args.output_root,
                    args.cohort,
                    rerender_report=args.rerender_report,
                    recompute_assessment=args.recompute_assessment,
                    case_inputs=read_patients(args.case_inputs, args.cohort) if args.case_inputs else (),
                )
            code = summary["exit_code"]
        print(json.dumps(summary, ensure_ascii=True, allow_nan=False, sort_keys=True))
        return code
    except MendelPrimeError as error:
        LOG.error("command=%s field=%s code=%s", args.command, error.field, error.code)
        return 1
    except KeyboardInterrupt:
        LOG.warning("command=%s interrupted", args.command)
        return 130
    except Exception as error:
        LOG.error("command=%s code=command_failed type=%s", args.command, type(error).__name__)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
