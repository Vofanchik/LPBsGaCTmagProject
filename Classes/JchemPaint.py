import subprocess
from .config import jchempaint_path


def runJCP():
    jar_file_path = jchempaint_path
    command = ["java", "-jar", jar_file_path]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
