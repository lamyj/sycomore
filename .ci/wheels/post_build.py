import glob
import site
import sys

subprocess.check_call([
    sys.executable, "-m", "pip", "install", "--user", *glob.glob("dist/*whl")])

sys.path.append(site.USER_SITE)
import sycomore
sycomore.gamma
