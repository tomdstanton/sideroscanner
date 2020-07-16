#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .tools.blast import run_blastp

def seed_blast(in_file, blastdb, eval, window, threads):
    hit_acc = []
    for q in run_blastp(in_file, blastdb, eval, threads):
        if len(q.hits) > 0:
            print(f'{len(q.hits)} hits found -> ', end="", flush=True)
            up = q.seq_len + window
            down = q.seq_len - window
            for h in q.hits:
                if down <= h.seq_len <= up:
                    hit_acc.append(h.blast_id)
            print(f"{len(hit_acc)} hits filtered between {down}-{up} aa")
    return hit_acc