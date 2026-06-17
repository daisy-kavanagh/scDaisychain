from scdaisychain.cli import main as pipeline_main


def test_pipeline_help_runs():
    assert pipeline_main([]) == 0
