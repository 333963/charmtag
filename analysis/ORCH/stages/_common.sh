# wspólne dla stages/*.sh (source przez pipeline.sh po załadowaniu etapu). Ścieżki jak w START_lxplus_*.command.
: "${CH:=$HOME/Documents/charm}"; W=$CH/charm_work; QL=$W/Q_lxplus; RQ=private/charm_q/Q_lxplus; EOSQ=/eos/user/k/klyczak/charm_q
queued() { rssh "condor_q -batch 2>/dev/null | grep -q '$1'"; }                         # 0 = batch w kolejce
rtest()  { rssh "test -f $1"; }
mdjson=(--prune-empty-dirs --include '*/' --include '*.md' --include '*.json' --include DONE --exclude 'cache/' --exclude 'models/' --exclude 'scores/' --exclude '*')   # tablica: bez globowania w cwd (BP 11-04)
