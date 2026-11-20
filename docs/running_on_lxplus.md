# Uruchamianie na lxplus (HTCondor, EOS, AFS)

Zasada: **Mac = kod, podpróbki (≤1,5 GB) i małe wyniki; lxplus = pełne pliki h5, treningi, pary, toys**. Wejścia i wyjścia na EOS, na AFS tylko kod i logi (kwota AFS ~2 GB — HTCondor odsyła pliki ze scratch joba do `initialdir`, więc duże pliki z joba muszą trafiać na EOS, a scratch ma być kasowany).

## Struktura paczki (`Q_lxplus/<etap>/`)
```
common.sh          # środowisko: LCG view / venv, PYTHONPATH (własna kopia src, nie editable), zmienne EOS_OUT, RUNS
<etap>_job.sh      # jeden job: mktemp -d w scratch, trap 'rm -rf $W' EXIT, wejścia z EOS, wyniki → EOS
<etap>.sub         # HTCondor submit (request_cpus/memory, +JobFlavour), should_transfer_files=NO
make_dag_<etap>.sh # generuje DAG (zależności: pairs → fit2d → toy×N → merge; syst; mix; compare)
submit_<etap>.sh   # check | submit | status | fetch
README.md
```
- `check`: `condor_q` (nie dublować biegu; główny DAG nie w kolejce), obecność wejść na EOS, kwota AFS <85 %, brak `DONE`.
- `submit`: `condor_submit_dag -batch-name charm_<etap>`; `status`: `condor_q -dag`, `DONE` na EOS; `fetch`: rsync wyników (`*.md`, `*.json`, bez h5/modeli) na Maca.
- venv (zuko, xgboost, iminuit) instalować **raz przed** DAG, nigdy `pip install` w jobie (zapis na AFS, wyścig węzłów).
- Toys dominują czas: DAG par z 50 toy-jobami ≈ 6–8 h; `NOTOYS=1` ≈ 3 h.

## Jeden dwuklik: `START_pipeline.command` (ORCH)
Orkiestrator na Macu (bash 3.2) sekwencjonuje istniejące `submit_*.sh` przez jedno połączenie ssh (ControlMaster):
- `pipeline.yaml`: etapy (`name`, `needs`, `max_hours`), np. `h2`, `pf_pairs`, `v`, `s5s8`, `s4`; `MAX_PARALLEL=2`.
- pętla co 10 min: `sync` (rsync paczki) → `check` → `submit` → `status` → po `DONE` `fetch`; stan w `state/`, raport `report.md`, log `pipeline.log`; powiadomienia macOS; `pipeline.stop` zatrzymuje; ponowny dwuklik wznawia.
- `pipeline.sh status|stop|reset <etap>`; `DRYRUN=1` ze stubem ssh do testów.

## 2FA na lxplus — co działa, co nie (test 04.11.2026)
- lxplus wymaga **drugiego czynnika przy każdym nowym połączeniu ssh** (PAM), niezależnie od biletu Kerberos. Klucze publiczne są odrzucane.
- Bilet Kerberos (25 h, odnawialny 5 dni, `kinit -R`) **nie omija 2FA**: GSSAPI daje `partial success`, potem żądanie OTP. Bilet z `kinit` (macOS Heimdal i MIT z Homebrew) nie niesie wskaźnika 2FA.
- Jedyne, co trwa godzinami, to **multipleksowanie ssh** (`ControlMaster auto`, `ControlPersist`): jeden OTP na sesję, kolejne `ssh`/`rsync` bez pytania. Sesja żyje tylko w Terminalu Maca; sandboxy (VM, chmura) ubijają procesy w tle, więc tam ControlMaster nie przeżyje.
- Zalecenie na dłużej: **acrontab na lxplus** (`check/submit/fetch` co 15 min po stronie CERN, wyniki na EOS) — wtedy Mac nie potrzebuje ssh poza pierwszą instalacją.

## Typowe awarie
| objaw | przyczyna | naprawa |
|---|---|---|
| joby HELD, DAG stoi | kwota AFS 100 % (pliki ze scratch odesłane do AFS) | przenieść logi/wyniki na EOS, `condor_release`; job kasuje własne kopie po `put_output` |
| `mux_client_request_session: Broken pipe` | ControlMaster zerwany | ponowne logowanie (OTP) i restart pipeline |
| PRE-skrypty DAG wiecznie DEFER | DAGMan na schedd nie widzi EOS / traci token | warunki PRE sprawdzać w jobie, nie w PRE; `DEFER` krótszy; katalog DAG poza AFS |
| xgboost nie importuje się | brak libomp / koło | fallback sklearn HistGB (inne liczby!) lub VENV_FIRST |
