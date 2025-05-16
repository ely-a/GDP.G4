from skyfield.api import load
import matplotlib.pyplot as plt

# Load ephemeris data and timescale
planets = load('de421.bsp')
ts = load.timescale()

# Define observation time
t = ts.utc(2033, 3, 10)

# Planet names and positions
planet_names = ['mercury', 'venus', 'earth', 'mars',
                'jupiter barycenter', 'saturn barycenter',
                'uranus barycenter', 'neptune barycenter']
positions = {}

for name in planet_names:
    planet = planets[name]
    pos = planet.at(t).ecliptic_position().au
    positions[name.title()] = pos

# Plot the solar system
fig, ax = plt.subplots(figsize=(8, 8))
for name, pos in positions.items():
    ax.plot(pos[0], pos[1], 'o', label=name.split()[0])
    ax.text(pos[0], pos[1], name.split()[0].capitalize(), fontsize=9, ha='right')

# Plot the Sun
ax.plot(0, 0, 'yo', label='Sun')

ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_title('Solar System Configuration - March 15, 2033')
ax.axis('equal')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()