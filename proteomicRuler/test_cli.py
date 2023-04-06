import os

from click.testing import CliRunner
from proteomicRuler.cli import main
from uuid import uuid4
class TestCli:
    def test_main(self):
        runner = CliRunner()
        filename = f"{uuid4()}.dat"
        result = runner.invoke(main, [
            "-i", r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For copy number test\combined\txt\proteinGroups.txt",
            "-o", filename,
            "-a", "Majority protein IDs",
            "-g"])
        os.remove(filename)
        assert result.exit_code == 0


