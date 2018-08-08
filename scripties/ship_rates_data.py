

import shutil
import os
import sys


package_name = sys.argv[1]
os.mkdir('../for_shipping/{0}'.format(package_name))

for sim_name in sys.argv[2:]:
    shutil.copy2('../runs/{0}/assoc_bootstrap.log'.format(sim_name),'../for_shipping/{0}/{1}.assoc_bootstrap.log'.format(package_name,sim_name))
    shutil.copy2('../runs/{0}/rates_bootstrap'.format(sim_name),'../for_shipping/{0}/{1}.rates_bootstrap'.format(package_name,sim_name))

