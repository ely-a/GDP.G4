from pykep.planet import jpl_lp

earth = jpl_lp('earth')
print(earth.name)         # Output: "Earth"
print(earth.name.lower()) # Output: "earth"
