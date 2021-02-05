#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sideroscanner.tools.blast import run_blastp

def seed_blast(in_file, blastdb, eval, window, threads):
    hit_acc = []
    for q in run_blastp(in_file, blastdb, eval, '10000', threads):
        if len(q.hits) > 0:
            print(f'{len(q.hits)} hits found -> ', end="", flush=True)
            up = q.seq_len + window
            down = q.seq_len - window
            for h in q.hits:
                if down <= h.seq_len <= up:
                    hit_acc.append(f'{h.blast_id}${h.seq_len}${"{:.1e}".format(h.hsps[0].evalue)}')
            if not hit_acc:
                print("\n[!] No hits for seed within size window, consider increasing window")
                return None
            else:
                print(f"{len(hit_acc)} hits filtered between {down}-{up} aa")
        else:
            print("\n[!] No blastp hits for seed, consider reducing E-value")
            return None
    return hit_acc
