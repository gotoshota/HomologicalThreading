import lammps_io as io
import homcloud.interface as hc
import numpy as np

lammpsdata = "../../../murashimaPolym/N2000/msd.N2000.1.data"
data = io.LammpsData(lammpsdata)

data.polyWrap()
coords = np.array(data.atoms.coords[0:1999])
pd = hc.PDList.from_alpha_filtration(coords)
pd1 = pd.dth_diagram(1)

print(pd1.births)

