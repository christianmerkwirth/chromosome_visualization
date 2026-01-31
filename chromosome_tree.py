import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
import numpy as np

# Configuration
CHROMOSOME_PAIRS = 3  # Number of autosomal pairs to simulate + Sex chromosome
CHROMOSOME_LENGTH = 100
CROSSOVER_RATE = 1.5  # Average crossovers per chromosome
Y_SCALE = 0.4 # Y chromosome is shorter

# Colors for Grandparents (Maternal GM, Maternal GF, Paternal GM, Paternal GF)
# Using distinctive colors: Red, Orange, Blue, Green
GP_COLORS = {
    'MGM': '#FF9999', # Maternal Grandmother (Reddish)
    'MGF': '#FFCC99', # Maternal Grandfather (Orangeish)
    'PGM': '#99CCFF', # Paternal Grandmother (Blueish)
    'PGF': '#99FF99'  # Paternal Grandfather (Greenish)
}

# Darker borders for contrast
GP_BORDERS = {
    'MGM': '#CC0000',
    'MGF': '#CC6600',
    'PGM': '#0000CC',
    'PGF': '#006600'
}

class Chromosome:
    def __init__(self, segments, is_y=False):
        # Segments is a list of (start, end, origin_id)
        self.segments = segments
        self.is_y = is_y

    def length(self):
        return CHROMOSOME_LENGTH * Y_SCALE if self.is_y else CHROMOSOME_LENGTH

def create_pure_chromosome(origin, is_y=False):
    return Chromosome([(0, CHROMOSOME_LENGTH, origin)], is_y)

def recombine(chr1, chr2):
    """
    Simulates meiosis with crossover between two homologous chromosomes.
    Returns ONE recombinant chromosome (haploid).
    """
    # If distinct sex chromosomes (XY), no recombination in this simple model (pseudoautosomal ignored)
    if chr1.is_y != chr2.is_y:
        # Return either X or Y with 50% prob
        return chr1 if random.random() < 0.5 else chr2

    # Standard Autosomal/XX Crossover
    length = CHROMOSOME_LENGTH
    # Determine number of crossovers (Poisson distribution usually, simplified here)
    num_crossovers = np.random.poisson(CROSSOVER_RATE)
    points = sorted([random.randint(0, length) for _ in range(num_crossovers)])
    points = [0] + points + [length]

    new_segments = []
    
    # Pick starting chromosome randomly
    current_source = chr1 if random.random() < 0.5 else chr2
    other_source = chr2 if current_source == chr1 else chr1
    
    for i in range(len(points) - 1):
        start, end = points[i], points[i+1]
        
        # Find segments in current_source that overlap with start, end
        # This is a simplification; we just take the color of the current source for this block
        # Assuming source is "pure" or already segmented.
        # To handle complex existing segments:
        
        # Helper to extract color map from source within range
        chunk_start = start
        while chunk_start < end:
            # Find which segment in source covers chunk_start
            src_seg = next(s for s in current_source.segments if s[0] <= chunk_start < s[1])
            seg_end = min(end, src_seg[1])
            new_segments.append((chunk_start, seg_end, src_seg[2]))
            chunk_start = seg_end
            
        # Swap strands
        current_source, other_source = other_source, current_source

    # Merge adjacent segments of same color
    merged = []
    if new_segments:
        curr = new_segments[0]
        for next_seg in new_segments[1:]:
            if next_seg[2] == curr[2] and next_seg[0] == curr[1]:
                curr = (curr[0], next_seg[1], curr[2])
            else:
                merged.append(curr)
                curr = next_seg
        merged.append(curr)
    
    return Chromosome(merged, is_y=chr1.is_y)

class Person:
    def __init__(self, name, gender, genome=None):
        self.name = name
        self.gender = gender # 'M' or 'F'
        self.genome = genome if genome else [] # List of pairs [ (chrA, chrB), ... ]

    def make_gamete(self):
        gamete = []
        for pair in self.genome:
            # For XY male, the last pair is XY.
            # If we are strictly following biology, X and Y don't recombine significantly for color tracking purposes here.
            # So we select either X or Y whole.
            # For autosomes, we recombine.
            if pair[0].is_y != pair[1].is_y:
                 # Sex chromosome pair in Male (X, Y)
                 # Randomly pick one
                 picked = pair[0] if random.random() < 0.5 else pair[1]
                 gamete.append(picked)
            else:
                # Homologous pair (Autosome or XX)
                gamete.append(recombine(pair[0], pair[1]))
        return gamete

def reproduce(mother, father, child_gender, child_name):
    mom_gamete = mother.make_gamete()
    dad_gamete = father.make_gamete()
    
    child_genome = []
    for i in range(len(mom_gamete)):
        m_chr = mom_gamete[i]
        f_chr = dad_gamete[i]
        
        # Ensure correct sex chromosome logic for final pair
        if i == len(mom_gamete) - 1: # Last pair is sex chromosomes
            # Mom always gives X. Dad gives X (daughter) or Y (son).
            # The gamete generation for Dad is random 50/50 X or Y.
            # We need to force it to match the requested child gender.
            
            # Find Dad's X and Y from his genome to know which color is which
            dad_sex_pair = father.genome[-1]
            dad_Y = dad_sex_pair[1] if dad_sex_pair[1].is_y else dad_sex_pair[0]
            dad_X = dad_sex_pair[0] if dad_sex_pair[1].is_y else dad_sex_pair[1]
            
            # Mom always gives X
            # (We use the one generated by recombination, which is fine, it's an X)
            
            if child_gender == 'M':
                # Force Dad's contribution to be the Y
                # We disregard the random dad_gamete selection for this specific chromosome
                f_chr = dad_Y 
            else:
                # Force Dad's contribution to be the X
                f_chr = dad_X 

        child_genome.append([m_chr, f_chr])
        
    return Person(child_name, child_gender, child_genome)


# --- Setup Generation 1 (Grandparents) ---
# Each has 3 pairs.
# MGM: XX, pure MGM color
# MGF: XY, pure MGF color
# PGM: XX, pure PGM color
# PGF: XY, pure PGF color

def make_founder(name, gender, color):
    genome = []
    for i in range(CHROMOSOME_PAIRS):
        is_sex_pair = (i == CHROMOSOME_PAIRS - 1)
        
        c1 = create_pure_chromosome(color, is_y=False) # Always X or Autosome
        
        if is_sex_pair and gender == 'M':
            c2 = create_pure_chromosome(color, is_y=True) # Y
        else:
            c2 = create_pure_chromosome(color, is_y=False) # X or Autosome
            
        genome.append([c1, c2])
    return Person(name, gender, genome)

MGM = make_founder("MGM", 'F', 'MGM')
MGF = make_founder("MGF", 'M', 'MGF')
PGM = make_founder("PGM", 'F', 'PGM')
PGF = make_founder("PGF", 'M', 'PGF')

# --- Setup Generation 2 (Parents) ---
# Mother: From MGM and MGF
# Father: From PGM and PGF

# To visualize "inheritance from grandparents", we create the parents by taking one pure copy from each grandparent
# (Simulating that parents are F1 of pure lines)

def make_parent(p1, p2, gender, name):
    # p1 is mom-side grandparent, p2 is dad-side grandparent
    # This function is manual to ensure they get distinct colored chromosomes without recombination "fuzz" yet,
    # or we can use the reproduce function. Let's use reproduce to be consistent, 
    # but since founders are pure, recombination doesn't change colors, just mixes "Red" with "Red".
    return reproduce(p1, p2, gender, name)

Mother = make_parent(MGM, MGF, 'F', "Mother")
Father = make_parent(PGM, PGF, 'M', "Father")

# --- Setup Generation 3 (Siblings) ---
# They will show recombination of colors.
Son = reproduce(Mother, Father, 'M', "Son")
Daughter = reproduce(Mother, Father, 'F', "Daughter")


# --- Plotting ---

fig, ax = plt.subplots(figsize=(12, 10))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Layout coordinates (approx percentage)
positions = {
    'MGM': (15, 85), 'MGF': (35, 85), 'PGM': (65, 85), 'PGF': (85, 85),
    'Mother': (25, 55), 'Father': (75, 55),
    'Daughter': (35, 20), 'Son': (65, 20)
}

def draw_person_genome(ax, person, center_x, center_y):
    # Draw label
    ax.text(center_x, center_y + 8, person.name, ha='center', fontsize=12, fontweight='bold')
    
    pair_spacing = 3.5
    chr_width = 1.5
    scale_y = 0.12 # Scale height of 100 unit chromosome to drawing units
    
    total_width = (CHROMOSOME_PAIRS * pair_spacing)
    start_x = center_x - (total_width / 2) + (pair_spacing/2)
    
    for pair_idx, pair in enumerate(person.genome):
        pair_x = start_x + (pair_idx * pair_spacing)
        
        # Draw each chromosome in the pair
        for ch_idx, chrom in enumerate(pair):
            offset_x = -0.5 if ch_idx == 0 else 0.5
            x = pair_x + offset_x - (chr_width/2)
            base_y = center_y - 5 # Shift down slightly
            
            # Draw segments
            for seg in chrom.segments:
                seg_len = seg[1] - seg[0]
                
                if chrom.is_y:
                    # Y is shorter, align to top (centromere approx)
                    pass
                
                # Rectangle(xy, width, height)
                # seg[0] is start (0-100), seg[1] is end
                # Invert Y so 0 is top
                seg_y = base_y + (CHROMOSOME_LENGTH - seg[1]) * scale_y 
                if chrom.is_y:
                    # Shift Y visual down to align tops if desired, or align bottoms?
                    # Usually aligned by centromere, but let's align tops for visual comparison
                    seg_y = base_y + (CHROMOSOME_LENGTH - seg[1]) * scale_y 
                
                rect_h = seg_len * scale_y
                
                color = GP_COLORS[seg[2]]
                border = GP_BORDERS[seg[2]]
                
                # Hatching/Texture for colorblind friendliness or extra distinction
                hatch = ''
                if seg[2] in ['MGF', 'PGF']:
                    hatch = '///'
                
                rect = patches.Rectangle((x, seg_y), chr_width, rect_h, 
                                         linewidth=0.5, edgecolor=border, facecolor=color, hatch=hatch)
                ax.add_patch(rect)
                
            # Outline
            total_h = chrom.length() * scale_y
            rect_y = base_y + (CHROMOSOME_LENGTH - chrom.length()) * scale_y if chrom.is_y else base_y
            
            # Highlight Y chromosome border for Male Lineage
            lw = 1.5 if chrom.is_y else 0.5
            ec = 'black'
            
            # Label Y
            if chrom.is_y:
                ax.text(x + chr_width/2, rect_y - 2, "Y", ha='center', fontsize=8, color='black')

# Draw lines connecting family tree
def draw_connections():
    # Grandparents to Parents
    # MGM+MGF -> Mother
    ax.plot([15, 25], [82, 65], color='black', alpha=0.3)
    ax.plot([35, 25], [82, 65], color='black', alpha=0.3)
    
    # PGM+PGF -> Father
    ax.plot([65, 75], [82, 65], color='black', alpha=0.3)
    ax.plot([85, 75], [82, 65], color='black', alpha=0.3)
    
    # Parents to Children
    # Center of parents
    parent_mid_x = 50
    parent_mid_y = 55
    
    # Line from Mom to Center
    ax.plot([25, 50], [48, 48], color='black', alpha=0.3)
    # Line from Dad to Center
    ax.plot([75, 50], [48, 48], color='black', alpha=0.3)
    
    # Line down from connection
    ax.plot([50, 50], [48, 35], color='black', alpha=0.3)
    
    # Split to kids
    ax.plot([50, 35], [35, 28], color='black', alpha=0.3)
    ax.plot([50, 65], [35, 28], color='black', alpha=0.3)


draw_connections()

for p_key, pos in positions.items():
    p_obj = eval(p_key)
    draw_person_genome(ax, p_obj, pos[0], pos[1])

# Legend
legend_elements = [
    patches.Patch(facecolor=GP_COLORS['MGM'], edgecolor=GP_BORDERS['MGM'], label='Maternal GM'),
    patches.Patch(facecolor=GP_COLORS['MGF'], edgecolor=GP_BORDERS['MGF'], hatch='///', label='Maternal GF'),
    patches.Patch(facecolor=GP_COLORS['PGM'], edgecolor=GP_BORDERS['PGM'], label='Paternal GM'),
    patches.Patch(facecolor=GP_COLORS['PGF'], edgecolor=GP_BORDERS['PGF'], hatch='///', label='Paternal GF (Y Source)'),
]
ax.legend(handles=legend_elements, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0))

plt.title("Random Recombination of Chromosomes in Siblings", fontsize=16)
plt.tight_layout()
plt.savefig("chromosome_inheritance.png", dpi=150)
print("Plot saved to chromosome_inheritance.png")
