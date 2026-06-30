import json
import logging
import os
import re

pjoin = os.path.join
logger = logging.getLogger(__name__)

# Match strings using one or more regular expressions
def _match_line(line, regex_lst, cfg_fpath):
    if len(regex_lst) == 0:
        return True
    for pat in regex_lst:
        try:
            m = re.search(r"{}".format(pat), line)
        except re.error as exc:
            raise ValueError(
                "Invalid MATCH regex {pat!r} while reading cfg {cfg_fpath!r}".format(
                    pat=pat,
                    cfg_fpath=cfg_fpath,
                )
            ) from exc
        if m is not None:
            return True
    return False


def regex_match(lst, regex_lst, cfg_fpath="<unknown cfg>"):
    # NOTE: For the regex_lst patterns, we use the raw string to generate the regular expression.
    #       This means that any regex special characters in the regex_lst should be properly
    #       escaped prior to calling this function.
    # NOTE: The input list is assumed to be a list of str objects and nothing else!
    matches = []
    for s in lst:
        if _match_line(s, regex_lst, cfg_fpath):
            matches.append(s)
    return matches


def load_json_file(fpath, cfg_fpath=None, cfg_line=None):
    if not os.path.exists(fpath):
        details = "JSON file not found: {fpath}".format(fpath=fpath)
        if cfg_fpath is not None:
            details += " while reading cfg {cfg_fpath}".format(cfg_fpath=cfg_fpath)
        if cfg_line is not None:
            details += " from cfg line {cfg_line!r}".format(cfg_line=cfg_line)
        raise FileNotFoundError(details)
    with open(fpath) as f:
        try:
            jsn = json.load(f)
        except json.JSONDecodeError as exc:
            details = "Invalid JSON in {fpath}".format(fpath=fpath)
            if cfg_fpath is not None:
                details += " while reading cfg {cfg_fpath}".format(cfg_fpath=cfg_fpath)
            if cfg_line is not None:
                details += " from cfg line {cfg_line!r}".format(cfg_line=cfg_line)
            raise ValueError(details) from exc
    for i, fn in enumerate(jsn.get('files', [])):
        fn = fn.replace("//", "/")
        jsn['files'][i] = fn
    return jsn


def read_cfg(fpath, match=[]):
    cfg_dir, fname = os.path.split(fpath)
    cfg_dir = cfg_dir or os.curdir
    if not os.path.exists(fpath):
        raise FileNotFoundError("Cfg file not found: {fpath}".format(fpath=fpath))
    if not os.path.exists(cfg_dir):
        raise FileNotFoundError("Cfg directory not found: {cfg_dir}".format(cfg_dir=cfg_dir))
    match = list(match)
    cfg = {
        "cfg_dir": cfg_dir,
        "src_xrd": "",
        "dst_xrd": "",
        "jsons": {},
    }
    sample_paths = {}
    blank_lines = 0
    comment_lines = 0
    metadata_lines = 0
    skipped_lines = 0
    with open(fpath) as f:
        for line_number, raw_line in enumerate(f, start=1):
            stripped_line = raw_line.strip()
            if not stripped_line:
                blank_lines += 1
                continue
            if stripped_line.startswith("#"):
                comment_lines += 1
                continue
            l = stripped_line.split("#", 1)[0].strip()
            if not len(l):
                comment_lines += 1
                continue
            if l.startswith("root:"):
                cfg['src_xrd'] = l
                metadata_lines += 1
            else:
                if len(regex_match([l], regex_lst=match, cfg_fpath=fpath)) == 0:
                    skipped_lines += 1
                    continue
                sample = os.path.basename(l)
                sample = sample.replace(".json", "")

                full_path = l if os.path.isabs(l) else os.path.normpath(pjoin(cfg['cfg_dir'], l))
                if sample in cfg['jsons']:
                    raise ValueError(
                        "Duplicate sample key {sample!r} while reading cfg {fpath!r}: "
                        "{first_path!r} and {second_path!r}".format(
                            sample=sample,
                            fpath=fpath,
                            first_path=sample_paths[sample],
                            second_path=full_path,
                        )
                    )
                jsn = load_json_file(full_path, cfg_fpath=fpath, cfg_line=l)
                cfg['jsons'][sample] = jsn
                sample_paths[sample] = full_path
                logger.debug(
                    "Selected cfg sample %s from %s line %s -> %s",
                    sample,
                    fpath,
                    line_number,
                    full_path,
                )
    logger.debug(
        "Read cfg %s: selected=%s blank=%s comments=%s metadata=%s skipped_by_match=%s match=%s",
        fpath,
        len(cfg['jsons']),
        blank_lines,
        comment_lines,
        metadata_lines,
        skipped_lines,
        match,
    )
    return cfg
