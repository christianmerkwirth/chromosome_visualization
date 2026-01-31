import requests
import json
import matplotlib.pyplot as plt
import numpy as np

def load_geojson(url):
    print(f"Downloading {url}...")
    try:
        r = requests.get(url)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"Failed to download: {e}")
        return None

def get_centroid_lon(geometry):
    # approximate centroid longitude from bounding box
    lons = []
    
    def extract_lons(coords):
        # coords can be [x, y] or list of lists
        # recursively find numbers
        if isinstance(coords[0], (float, int)):
            lons.append(coords[0])
        else:
            for c in coords:
                extract_lons(c)
                
    extract_lons(geometry['coordinates'])
    
    if not lons:
        return 0
    return (min(lons) + max(lons)) / 2

def plot_geometry(ax, geometry, color='k', linewidth=0.5, shift_x=0):
    geo_type = geometry['type']
    coords = geometry['coordinates']
    
    def plot_poly(poly_coords):
        # poly_coords is list of rings. First is exterior.
        exterior = poly_coords[0]
        xs = [pt[0] + shift_x for pt in exterior]
        ys = [pt[1] for pt in exterior]
        ax.plot(xs, ys, color=color, linewidth=linewidth)
        
    if geo_type == 'Polygon':
        plot_poly(coords)
    elif geo_type == 'MultiPolygon':
        for poly in coords:
            plot_poly(poly)

# URLs
WORLD_URL = "https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json"
US_STATES_URL = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"

# Load Data
world_data = load_geojson(WORLD_URL)
states_data = load_geojson(US_STATES_URL)

if not world_data or not states_data:
    print("Error downloading data. Exiting.")
    exit(1)

# Find Algeria and Texas
algeria_geom = None
texas_geom = None

for feature in world_data['features']:
    if feature['properties'].get('name') == 'Algeria':
        algeria_geom = feature['geometry']
        break
        
for feature in states_data['features']:
    # Check simplified vs regular properties
    props = feature['properties']
    if props.get('name') == 'Texas' or props.get('NAME') == 'Texas':
        texas_geom = feature['geometry']
        break

if not algeria_geom:
    print("Could not find Algeria.")
if not texas_geom:
    print("Could not find Texas.")
    
if not algeria_geom or not texas_geom:
    exit(1)

# Calculate Shift
algeria_lon = get_centroid_lon(algeria_geom)
texas_lon = get_centroid_lon(texas_geom)
base_shift = algeria_lon - texas_lon

# Extra shift: 200 km to the east
# Approx conversion at ~32 deg latitude: 1 deg lon ~= 94.4 km
# 200 km / 94.4 km/deg ~= 2.12 deg
extra_shift_deg = 200 / 94.4
shift = base_shift + extra_shift_deg

print(f"Base Shift: {base_shift:.2f} degrees")
print(f"Extra Shift (200km): {extra_shift_deg:.2f} degrees")
print(f"Total Shift: {shift:.2f} degrees")

# Plot
fig, ax = plt.subplots(figsize=(20, 12))
ax.set_aspect('equal')
ax.axis('off')

# Plot World
print("Plotting World...")
for feature in world_data['features']:
    # Skip Antarctica for cleaner view? Optional.
    if feature['properties'].get('name') == 'Antarctica':
        continue
    plot_geometry(ax, feature['geometry'], color='black', linewidth=0.5)

# Plot Shifted Texas
print("Plotting Shifted Texas...")
plot_geometry(ax, texas_geom, color='red', linewidth=1.5, shift_x=shift)

# Title
plt.title("Texas Transposed to Algeria's Longitude (Latitude Preserved)", fontsize=16)

# Zoom into Northern Africa / Mediterranean
ax.set_xlim(-20, 40)
ax.set_ylim(15, 50)

# Save
output_file = "texas_over_algeria.png"
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"Saved map to {output_file}")
