import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
import numpy as np

# Configuration
# 9 Autosomes + 1 Sex Pair = 10 Pairs
CHROMOSOME_PAIRS = 10 
# Relative lengths for first 9 + X/Y
CHROMOSOME_LENGTHS = {
    1: 249, 2: 242, 3: 198, 4: 190, 5: 181, 6: 171, 7: 159, 8: 145, 9: 138, 
    'X': 156, 'Y': 57
}
CROSSOVER_RATE_FACTOR = 0.012

# Colors for Grandparents
GP_COLORS = {
    'MGM': '#FF9999', # Pink
    'MGF': '#FFCC99', # Orange
    'PGM': '#99CCFF', # Blue
    'PGF': '#99FF99'  # Green
}

GP_BORDERS = {
    'MGM': '#CC0000',
    'MGF': '#CC6600',
    'PGM': '#0000CC',
    'PGF': '#006600'
}

# Textures for the two homologous chromosomes of a grandparent
# 0 = "Left" chromosome (Solid)
# 1 = "Right" chromosome (Hatched)
HATCH_PATTERNS = {
    0: '',      # Solid
    1: '....'   # Dotted/Stippled
}

class Chromosome:
    def __init__(self, segments, name, length):
        # Segments: list of (start, end, origin_data)
        # origin_data is (gp_id, chromatid_id)
        self.segments = segments
        self.name = name 
        self.total_length = length

    def length(self):
        return self.total_length

def create_pure_chromosome(name, gp_id, chromatid_id):
    # gp_id: 'MGM', etc.
    # chromatid_id: 0 or 1 (to distinguish the pair)
    length = CHROMOSOME_LENGTHS[name]
    return Chromosome([(0, length, (gp_id, chromatid_id))], name, length)

def recombine(chr1, chr2):
    # XY logic: No recombination, return one or the other
    if chr1.name != chr2.name:
         return chr1 if random.random() < 0.5 else chr2

    length = chr1.length()
    expected_crossovers = length * CROSSOVER_RATE_FACTOR
    num_crossovers = np.random.poisson(expected_crossovers)
    
    points = sorted([random.uniform(0, length) for _ in range(num_crossovers)])
    points = [0] + points + [length]

    new_segments = []
    
    current_source = chr1 if random.random() < 0.5 else chr2
    other_source = chr2 if current_source == chr1 else chr1
    
    for i in range(len(points) - 1):
        start, end = points[i], points[i+1]
        
        # Extract slices from current_source
        chunk_start = start
        while chunk_start < end:
            src_seg = next(s for s in current_source.segments if s[0] <= chunk_start < s[1])
            seg_end = min(end, src_seg[1])
            # Keep the origin data (gp_id, chromatid_id)
            new_segments.append((chunk_start, seg_end, src_seg[2]))
            chunk_start = seg_end
            
        current_source, other_source = other_source, current_source

    return Chromosome(new_segments, chr1.name, length)

class Person:
    def __init__(self, name, gender, genome=None):
        self.name = name
        self.gender = gender 
        self.genome = genome if genome else [] 

    def make_gamete(self):
        gamete = []
        for pair in self.genome:
            if pair[0].name != pair[1].name:
                 # XY pair
                 picked = pair[0] if random.random() < 0.5 else pair[1]
                 gamete.append(picked)
            else:
                # Homologous pair
                gamete.append(recombine(pair[0], pair[1]))
        return gamete

def reproduce(mother, father, child_gender, child_name):
    mom_gamete = mother.make_gamete()
    dad_gamete = father.make_gamete()
    
    child_genome = []
    
    # We are simulating 9 autosomes (indices 0-8) + 1 sex pair (index 9)
    for i in range(len(mom_gamete)):
        m_chr = mom_gamete[i]
        
        # Sex chromosome logic
        if i == CHROMOSOME_PAIRS - 1: # Last one
            dad_sex_pair = father.genome[-1]
            dad_X = dad_sex_pair[0] if dad_sex_pair[0].name == 'X' else dad_sex_pair[1]
            dad_Y = dad_sex_pair[0] if dad_sex_pair[0].name == 'Y' else dad_sex_pair[1]
            
            if child_gender == 'M':
                f_chr = dad_Y 
            else:
                f_chr = dad_X 
        else:
            f_chr = dad_gamete[i]
            
        child_genome.append([m_chr, f_chr])
        
    return Person(child_name, child_gender, child_genome)

# --- Setup Generation 1 (Founders) ---
def make_founder(name, gender, gp_id):
    genome = []
    # Autosomes 1-9
    for i in range(1, 10):
        # Create pair with distinct textures (0 and 1)
        c1 = create_pure_chromosome(i, gp_id, 0)
        c2 = create_pure_chromosome(i, gp_id, 1)
        genome.append([c1, c2])
        
    # Sex Chromosomes
    if gender == 'M':
        cx = create_pure_chromosome('X', gp_id, 0)
        cy = create_pure_chromosome('Y', gp_id, 1) # Y gets texture 1
        genome.append([cx, cy])
    else:
        c1 = create_pure_chromosome('X', gp_id, 0)
        c2 = create_pure_chromosome('X', gp_id, 1)
        genome.append([c1, c2])
        
    return Person(name, gender, genome)

MGM = make_founder("MGM (Maryna)", 'F', 'MGM')
MGF = make_founder("MGF (Valery)", 'M', 'MGF')
PGM = make_founder("PGM (Marlene)", 'F', 'PGM')
PGF = make_founder("PGF (Helmuth)", 'M', 'PGF')

# --- Setup Generation 2 (Parents) ---
Mother = reproduce(MGM, MGF, 'F', "Mother (Svitlana)")
Father = reproduce(PGM, PGF, 'M', "Father (Christian)")

# --- Setup Generation 3 (Siblings) ---
Son = reproduce(Mother, Father, 'M', "Son (Hugo)")
Daughter = reproduce(Mother, Father, 'F', "Daughter (Viki)")

# --- Plotting ---
fig, ax = plt.subplots(figsize=(16, 12))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

def draw_person_genome(ax, person, center_x, center_y, width_limit):
    ax.text(center_x, center_y + 8, person.name, ha='center', fontsize=12, fontweight='bold')
    
    num_pairs = CHROMOSOME_PAIRS
    pair_spacing = width_limit / num_pairs 
    chr_width = pair_spacing * 0.4
    scale_y = 0.05
    
    start_x = center_x - (width_limit / 2) + (pair_spacing/2)
    
    for pair_idx, pair in enumerate(person.genome):
        pair_x = start_x + (pair_idx * pair_spacing)
        
        # Label
        chr_name = str(pair_idx + 1) if pair_idx < 9 else "X/Y"
        ax.text(pair_x, center_y - 12, chr_name, ha='center', fontsize=8, color='#333')
        
        for ch_idx, chrom in enumerate(pair):
            offset_x = - (chr_width * 0.55) if ch_idx == 0 else (chr_width * 0.55)
            x = pair_x + offset_x - (chr_width/2)
            
            base_y = center_y 
            
            for seg in chrom.segments:
                seg_start, seg_end, (gp_id, chromatid_id) = seg
                
                rect_y = base_y - (seg_end * scale_y)
                rect_h = (seg_end - seg_start) * scale_y
                
                color = GP_COLORS[gp_id]
                border = GP_BORDERS[gp_id]
                hatch = HATCH_PATTERNS[chromatid_id]
                
                # Make hatch finer/denser if needed
                # For matplotlib patches, hatch density is controlled by repeating char, e.g. '||||'
                
                rect = patches.Rectangle((x, rect_y), chr_width, rect_h, 
                                         linewidth=0.5, edgecolor=border, facecolor=color, hatch=hatch)
                ax.add_patch(rect)
                
            if chrom.name == 'Y':
                ax.text(x + chr_width/2, base_y - (chrom.length()*scale_y) - 3, "Y", ha='center', fontsize=7, fontweight='bold')

# Layout
row1_y = 85
row2_y = 55
row3_y = 20

gp_width = 20
p_width = 35
c_width = 35

draw_person_genome(ax, MGM, 15, row1_y, gp_width)
draw_person_genome(ax, MGF, 35, row1_y, gp_width)
draw_person_genome(ax, PGM, 65, row1_y, gp_width)
draw_person_genome(ax, PGF, 85, row1_y, gp_width)

draw_person_genome(ax, Mother, 25, row2_y, p_width)
draw_person_genome(ax, Father, 75, row2_y, p_width)

draw_person_genome(ax, Daughter, 25, row3_y, c_width)
draw_person_genome(ax, Son, 75, row3_y, c_width)

# Connections
def connect(x1, y1, x2, y2):
    ax.plot([x1, x2], [y1, y2], color='black', alpha=0.2, linewidth=1)

connect(15, row1_y-10, 25, row2_y+10)
connect(35, row1_y-10, 25, row2_y+10)
connect(65, row1_y-10, 75, row2_y+10)
connect(85, row1_y-10, 75, row2_y+10)

connect(25, row2_y-10, 25, row3_y+10)
connect(75, row2_y-10, 25, row3_y+10)
connect(25, row2_y-10, 75, row3_y+10)
connect(75, row2_y-10, 75, row3_y+10)

# Legend
legend_elements = [
    patches.Patch(facecolor='white', edgecolor='black', label='Grandparents:'),
    patches.Patch(facecolor=GP_COLORS['MGM'], edgecolor=GP_BORDERS['MGM'], label='Maternal GM'),
    patches.Patch(facecolor=GP_COLORS['MGF'], edgecolor=GP_BORDERS['MGF'], label='Maternal GF'),
    patches.Patch(facecolor=GP_COLORS['PGM'], edgecolor=GP_BORDERS['PGM'], label='Paternal GM'),
    patches.Patch(facecolor=GP_COLORS['PGF'], edgecolor=GP_BORDERS['PGF'], label='Paternal GF'),
    patches.Patch(facecolor='white', edgecolor='black', label=' | '), # Spacer
    patches.Patch(facecolor='lightgrey', edgecolor='black', label='Textures (Homologs):'),
    patches.Patch(facecolor='white', edgecolor='black', label='Chr A (Solid)'),
    patches.Patch(facecolor='white', edgecolor='black', hatch='....', label='Chr B (Stippled)'),
]
ax.legend(handles=legend_elements, loc='lower center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 0))

plt.suptitle("Chromosome Inheritance: 10 Pairs with Recombination & Homolog Tracking", fontsize=18)
plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.savefig("chromosome_inheritance_10pairs.png", dpi=150)
print("Saved to chromosome_inheritance_10pairs.png")
