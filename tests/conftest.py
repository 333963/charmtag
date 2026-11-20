"""Testy wymagające pakietu `charmbkg` (projekt B, prywatne repo charmbkg) są pomijane, gdy pakietu nie ma.

`charmtag` zależy od `charmbkg` tylko w ścieżkach projektu B (aliasy `charmtag.bkggen`, `charmtag.hyp`, … oraz źródła tła
treningowego inne niż sideband). GitHub Actions nie ma dostępu do prywatnego repo, więc tam działa wyłącznie część A;
lokalnie z zainstalowanym `charmbkg` uruchamia się wszystko. Pominięcie jest jawne: nagłówek i podsumowanie pytest.
"""
import importlib.util

HAVE_CHARMBKG = importlib.util.find_spec("charmbkg") is not None

# pliki, które importują moduły projektu B już przy kolekcji
B_TEST_FILES = ["test_aa_merge.py", "test_b4.py", "test_bkggen.py", "test_candrot.py", "test_dstar_tag.py", "test_feta.py", "test_hyp.py",
                "test_s5_bkg.py", "test_s8_salad.py", "test_salad_v1.py"]

collect_ignore = [] if HAVE_CHARMBKG else list(B_TEST_FILES)


def _msg():
    return (f"charmbkg: BRAK — pominięto {len(B_TEST_FILES)} plików testów projektu B ({', '.join(B_TEST_FILES)}); "
            "pojedyncze testy A wymagające charmbkg oznaczone jako skipped")


def pytest_report_header(config):
    return "charmbkg: zainstalowany — pełny zestaw testów" if HAVE_CHARMBKG else _msg()


def pytest_terminal_summary(terminalreporter):
    if not HAVE_CHARMBKG:
        terminalreporter.write_line(_msg(), yellow=True)
