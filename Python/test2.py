from pykep.planet import jpl_lp

jupiter = jpl_lp('jupiter')

# Compute the minimum flyby distance in terms of planetary radii
min_flyby_radii = jupiter.safe_radius / jupiter.radius

print(f"Minimum flyby distance: {min_flyby_radii:.4f} Jupiter radii")
