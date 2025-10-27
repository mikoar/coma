from __future__ import annotations
import os, socket, math
from typing import TextIO, List
from src.alignment.alignment_results import SmapResultRow

class SmapReader:
    @staticmethod
    def format_field(v):
        if v is None: return "NA"
        if isinstance(v, float) and math.isnan(v): return "NA"
        return str(v)

    @staticmethod
    def arg_to_string(v):
        if v is None: return "None"
        if isinstance(v, (list, tuple)): return "[" + ", ".join(str(x) for x in v) + "]"
        return str(v)

    def write_headers(self, f: TextIO, ref_path: str, qry_path: str, args, version="0.9"):
        f.writelines(line + "\n" for line in [
            f"# hostname={socket.gethostname()}",
            "# coma " + " ".join([f"--{k} {self.arg_to_string(v)}" for k, v in vars(args).items()]),
            f"# SMAP File Version:\t{version}",
            f"# Reference Maps From:\t{os.path.abspath(ref_path) if ref_path else ''}",
            f"# Query Maps From:\t{os.path.abspath(qry_path) if qry_path else ''}",
        ])

    def writeSmapsFromRows(self, outputFile, smap_rows: List["SmapResultRow"], args):
        types = SmapResultRow.COLUMN_TYPES
        names = SmapResultRow.COLUMN_NAMES

        rows = [r.to_smap_dict() for r in smap_rows]

        for i, d in enumerate(rows, start=1):
            d["SmapEntryID"] = i
            if d.get("LinkID") is None:
                d["LinkID"] = -1

        for i, r in enumerate(smap_rows):
            peer_idx = getattr(r, "_link_peer_idx", None)
            if isinstance(peer_idx, int) and 0 <= peer_idx < len(rows):
                rows[i]["LinkID"] = rows[peer_idx]["SmapEntryID"]

        with open(outputFile.name if hasattr(outputFile, "name") else outputFile, "w", encoding="utf-8") as f:
            f.write("# SMAP File Version:\t0.9\n")
            f.write("#h " + "\t".join(names) + "\n")
            f.write("#f " + "\t".join(types) + "\n")

            for d in rows:
                out = []
                for col, typ in zip(names, types):
                    val = d.get(col, None)
                    if val is None:
                        out.append("NA")
                    else:
                        if typ == "int":
                            out.append(str(int(val)))
                        elif typ == "float":
                            out.append(str(float(val)))
                        else:
                            out.append(str(val))
                f.write("\t".join(out) + "\n")
