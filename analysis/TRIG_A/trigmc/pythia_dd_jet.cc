// pythia_dd_jet.cc — TRIG_A/trigmc (projekt A, 17.11.2026). Pochodna Q_lxplus/syst/pythia_dd_syst.cc (generator-level Pythia 8, A3),
// rozszerzona o wielkości potrzebne do EMULACJI TRIGGERA L1_J12 na poziomie cząstek i testu zamknięcia wag w2(n_trk, lead_pt):
//   * okno kalorymetryczne 0,8×0,8 (4×4 wieże 0,2×0,2, |eta|<3,2, przesuwne) — maksimum sumy E_T:
//       etem  = skala EM (gamma, e: 1,0; hadrony: HADRESP; naładowane pT<0,4 GeV i miony pominięte),
//       etall = wszystkie widzialne cząstki (bez neutrin i mionów), odpowiedź 1,0;
//   * dżety naładowane anty-kT R=0,4 (naładowane pT>0,5, |eta|<2,5) — tylko w zdarzeniach z D;
//   * histogram 3D WSZYSTKICH zdarzeń (n_ch, lead_pt, E_T okna) — mianownik wag typu w2.
// Użycie: pythia_dd_jet <nEvents> <seed> <out.tsv> [tune=A3|Monash] [eCM=13600] [procs=nd]
// Wyjście TSV (tylko zdarzenia z >=1 D o pT>2, |eta|<3):
//   E <ev> <proc> <nch> <lead> <pthat> <etem> <etall> <weta> <wphi> <cjlead>
//   D <ev> <pdg> <pt> <eta> <phi> <fromB> <mode> <cjpt>      (cjpt: pT dżetu naładowanego o osi w dR<0,4 od D; 0 = brak)
// + <out>.json (liczniki, sigma, ustawienia) + <out>.h3em / <out>.h3all („klucz liczba”, klucz = (nb*NL + lb)*NE + eb).
#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace Pythia8;

static const double PTSTORE = 2.0, ETASTORE = 3.0, HADRESP = 0.65;
static const int NN = 200, NL = 60, NE = 80;          // n_ch 0..199 (krok 1), lead 0..30 (0,5 GeV), E_T 0..40 (0,5 GeV); ostatni = przepełnienie
static const int NETA = 32, NPHI = 32;                // wieże 0,2 x 0,196

static bool isB(int id) { int a = std::abs(id); return (a / 100) % 10 == 5 || (a / 1000) % 10 == 5; }
static void finals(const Event& ev, int i, std::vector<int>& out, int depth = 0) {
  if (depth > 20) return;
  if (ev[i].isFinal()) { out.push_back(ev[i].idAbs()); return; }
  for (int k : ev[i].daughterList()) finals(ev, k, out, depth + 1);
}
static int modeOf(const Event& ev, int i) {
  const Particle& p = ev[i]; int a = p.idAbs();
  std::vector<int> fs; finals(ev, i, fs); std::sort(fs.begin(), fs.end());
  if (a == 421) return (fs == std::vector<int>{211, 321}) ? 1 : 0;
  if (a == 411) return (fs == std::vector<int>{211, 211, 321}) ? 1 : 0;
  if (a == 413) {
    std::vector<int> d; for (int k : p.daughterList()) d.push_back(ev[k].idAbs()); std::sort(d.begin(), d.end());
    if (d != std::vector<int>{211, 421}) return 0;
    return (fs == std::vector<int>{211, 211, 321}) ? 1 : 0;
  }
  return 0;
}
// maksimum sumy w oknie 4x4 (phi cykliczne); zwraca E_T i środek okna
static double maxWindow(const double g[NETA][NPHI], double& weta, double& wphi) {
  double best = 0.; int bi = 0, bj = 0;
  for (int i = 0; i + 4 <= NETA; ++i)
    for (int j = 0; j < NPHI; ++j) {
      double s = 0.;
      for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) s += g[i + a][(j + b) % NPHI];
      if (s > best) { best = s; bi = i; bj = j; }
    }
  weta = -3.2 + 0.2 * (bi + 2);
  wphi = -M_PI + (2. * M_PI / NPHI) * (bj + 2); if (wphi > M_PI) wphi -= 2. * M_PI;
  return best;
}
static void dumpH(const std::unordered_map<int, long>& h, const std::string& fn) {
  std::ofstream f(fn); for (auto& kv : h) f << kv.first << " " << kv.second << "\n";
}

int main(int argc, char** argv) {
  if (argc < 4) { fprintf(stderr, "uzycie: %s nEvents seed out.tsv [A3|Monash] [eCM] [nd,sd,dd,cd]\n", argv[0]); return 1; }
  long nEv = std::atol(argv[1]); int seed = std::atoi(argv[2]); std::string out = argv[3];
  std::string tune = argc > 4 ? argv[4] : "A3";
  double eCM = argc > 5 ? std::atof(argv[5]) : 13600.;
  std::string procs = argc > 6 ? argv[6] : "nd";

  Pythia pythia("", false);
  std::vector<std::string> S = {
      "Beams:idA = 2212", "Beams:idB = 2212", "Beams:eCM = " + std::to_string(eCM),
      "Random:setSeed = on", "Random:seed = " + std::to_string(seed % 900000000 + 1),
      "Next:numberCount = 0", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0",
      "Init:showChangedSettings = on", "Init:showChangedParticleData = off",
      std::string("SoftQCD:nonDiffractive = ") + (procs.find("nd") != std::string::npos ? "on" : "off"),
      std::string("SoftQCD:singleDiffractive = ") + (procs.find("sd") != std::string::npos ? "on" : "off"),
      std::string("SoftQCD:doubleDiffractive = ") + (procs.find("dd") != std::string::npos ? "on" : "off"),
      std::string("SoftQCD:centralDiffractive = ") + (procs.find("cd") != std::string::npos ? "on" : "off"),
      "ParticleDecays:limitTau0 = on", "ParticleDecays:tau0Max = 10.",
  };
  if (tune == "A3") {   // identycznie jak truth/pythia_dd.cc (baza normalizacji projektu A)
    std::vector<std::string> A3 = {
        "Tune:ee = 7", "Tune:pp = 14", "PDF:pSet = 13",
        "MultipartonInteractions:bProfile = 2", "MultipartonInteractions:pT0Ref = 2.45",
        "MultipartonInteractions:ecmPow = 0.21", "MultipartonInteractions:coreRadius = 0.55",
        "MultipartonInteractions:coreFraction = 0.9",
        "SigmaDiffractive:PomFlux = 4", "SigmaDiffractive:PomFluxEpsilon = 0.07",
        "SigmaDiffractive:PomFluxAlphaPrime = 0.25", "ColourReconnection:range = 1.8"};
    S.insert(S.end(), A3.begin(), A3.end());
  } else if (tune != "Monash") { fprintf(stderr, "nieznany tune %s\n", tune.c_str()); return 2; }
  if (procs.find("cd") != std::string::npos) S.push_back("SigmaTotal:zeroAXB = off");
  if (const char* ex = std::getenv("PYTHIA_EXTRA")) {
    std::stringstream ss(ex); std::string item;
    while (std::getline(ss, item, ';')) if (item.find_first_not_of(" ") != std::string::npos) S.push_back(item);
  }
  for (auto& s : S) if (!pythia.readString(s)) { fprintf(stderr, "BLAD ustawienia: %s\n", s.c_str()); return 3; }
  if (!pythia.init()) { fprintf(stderr, "BLAD init\n"); return 4; }

  SlowJet sj(-1, 0.4, 2.0, 2.5, 1, 1);            // anty-kT R=0.4, pT_jet>2, wszystkie czastki zdarzenia pomocniczego, masa pionu
  Event chg; chg.init("naladowane", &pythia.particleData);

  std::ofstream f(out);
  f << "# pythia_dd_jet seed=" << seed << " tune=" << tune << " eCM=" << eCM << " procs=" << procs << " hadresp=" << HADRESP << "\n";
  std::map<int, long> nProc;
  std::unordered_map<int, long> hem, hall;
  long nOk = 0, nFail = 0, nWithD = 0;
  static double gem[NETA][NPHI], gall[NETA][NPHI];
  for (long iEv = 0; iEv < nEv; ++iEv) {
    if (!pythia.next()) { if (++nFail > 1000) break; continue; }
    ++nOk;
    const Event& ev = pythia.event;
    int code = pythia.info.code(); nProc[code]++;
    int nch = 0; double lead = 0.;
    std::vector<int> ds;
    for (int i = 0; i < NETA; ++i) for (int j = 0; j < NPHI; ++j) { gem[i][j] = 0.; gall[i][j] = 0.; }
    for (int i = 0; i < ev.size(); ++i) {
      const Particle& p = ev[i];
      int a = p.idAbs();
      if (p.isFinal()) {
        double eta = p.eta(), pt = p.pT();
        if (p.isCharged() && pt > 0.5 && std::abs(eta) < 2.5) { ++nch; lead = std::max(lead, pt); }
        if (std::abs(eta) < 3.2 && a != 12 && a != 14 && a != 16 && a != 13) {
          int ie = std::min(NETA - 1, std::max(0, (int)std::floor((eta + 3.2) / 0.2)));
          int ip = std::min(NPHI - 1, std::max(0, (int)std::floor((p.phi() + M_PI) / (2. * M_PI / NPHI))));
          double et = p.eT();
          gall[ie][ip] += et;
          if (a == 22 || a == 11) gem[ie][ip] += et;
          else if (!(p.isCharged() && pt < 0.4)) gem[ie][ip] += HADRESP * et;
        }
      }
      if ((a == 421 || a == 411 || a == 413) && p.iBotCopyId() == i && p.pT() > PTSTORE && std::abs(p.eta()) < ETASTORE) ds.push_back(i);
    }
    double we, wp, we2, wp2;
    double etem = maxWindow(gem, we, wp), etall = maxWindow(gall, we2, wp2);
    int nb = std::min(nch, NN - 1), lb = std::min((int)(lead / 0.5), NL - 1);
    hem[(nb * NL + lb) * NE + std::min((int)(etem / 0.5), NE - 1)]++;
    hall[(nb * NL + lb) * NE + std::min((int)(etall / 0.5), NE - 1)]++;
    if (ds.empty()) continue;
    ++nWithD;
    chg.reset();
    for (int i = 0; i < ev.size(); ++i) {
      const Particle& p = ev[i];
      if (p.isFinal() && p.isCharged() && p.pT() > 0.5 && std::abs(p.eta()) < 2.5) chg.append(p);
    }
    sj.analyze(chg);
    double cjlead = sj.sizeJet() > 0 ? sj.pT(0) : 0.;
    char buf[320];
    snprintf(buf, sizeof buf, "E\t%ld\t%d\t%d\t%.4f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.3f\n", iEv, code, nch, lead, pythia.info.pTHat(), etem, etall, we, wp, cjlead);
    f << buf;
    for (int i : ds) {
      const Particle& p = ev[i];
      bool fromB = false;
      std::vector<int> stack = p.motherList(); std::vector<char> seen(ev.size(), 0); int guard = 0;
      while (!stack.empty() && guard++ < 2000) {
        int m = stack.back(); stack.pop_back();
        if (m <= 0 || seen[m]) continue;
        seen[m] = 1;
        if (isB(ev[m].id())) { fromB = true; break; }
        for (int mm : ev[m].motherList()) stack.push_back(mm);
      }
      double cj = 0.;
      for (int j = 0; j < sj.sizeJet(); ++j) {
        double dphi = std::abs(sj.phi(j) - p.phi()); if (dphi > M_PI) dphi = 2. * M_PI - dphi;
        double dr = std::sqrt(dphi * dphi + (sj.y(j) - p.y()) * (sj.y(j) - p.y()));
        if (dr < 0.4) { cj = sj.pT(j); break; }      // dzety posortowane malejaco w pT
      }
      snprintf(buf, sizeof buf, "D\t%ld\t%d\t%.5f\t%.5f\t%.5f\t%d\t%d\t%.3f\n", iEv, p.id(), p.pT(), p.eta(), p.phi(), (int)fromB, modeOf(ev, i), cj);
      f << buf;
    }
  }
  f.close();
  dumpH(hem, out + ".h3em"); dumpH(hall, out + ".h3all");
  std::ofstream js(out + ".json");
  js << "{\n \"seed\": " << seed << ", \"tune\": \"" << tune << "\", \"eCM\": " << eCM << ", \"procs\": \"" << procs << "\",\n";
  js << " \"n_requested\": " << nEv << ", \"n_ok\": " << nOk << ", \"n_fail\": " << nFail << ", \"n_with_D\": " << nWithD << ",\n";
  js << " \"hadresp\": " << HADRESP << ", \"NN\": " << NN << ", \"NL\": " << NL << ", \"NE\": " << NE << ",\n";
  js << " \"sigma_gen_mb_total\": " << pythia.info.sigmaGen() << ", \"sigma_err_mb_total\": " << pythia.info.sigmaErr() << ",\n";
  js << " \"proc\": {";
  bool first = true;
  for (auto& kv : nProc) { js << (first ? "" : ", ") << "\"" << kv.first << "\": {\"n\": " << kv.second << ", \"sigma_mb\": " << pythia.info.sigmaGen(kv.first) << "}"; first = false; }
  js << "},\n \"pythia_version\": " << PYTHIA_VERSION << "\n}\n";
  js.close();
  pythia.stat();
  return 0;
}
