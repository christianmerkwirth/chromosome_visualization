import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
import numpy as np

# Configuration
# Approx relative lengths of human chromosomes 1-22. 
# Data based on base pairs (millions).
CHROMOSOME_LENGTHS = {
    1: 249, 2: 242, 3: 198, 4: 190, 5: 181, 6: 171, 7: 159, 8: 145, 9: 138, 10: 133,
    11: 135, 12: 133, 13: 114, 14: 107, 15: 101, 16: 90, 17: 83, 18: 80, 19: 58, 20: 64,
    21: 46, 22: 50, 'X': 156, 'Y': 57
}

CHROMOSOME_PAIRS = 23
CROSSOVER_RATE_FACTOR = 0.012 # Crossovers per unit length. Long chroms have more.

# Colors for Grandparents (Maternal GM, Maternal GF, Paternal GM, Paternal GF)
GP_COLORS = {
    'MGM': '#FF9999', # Maternal Grandmother (Reddish)
    'MGF': '#FFCC99', # Maternal Grandfather (Orangeish)
    'PGM': '#99CCFF', # Paternal Grandmother (Blueish)
    'PGF': '#99FF99'  # Paternal Grandfather (Greenish)
}

GP_BORDERS = {
    'MGM': '#CC0000',
    'MGF': '#CC6600',
    'PGM': '#0000CC',
    'PGF': '#006600'
}

class Chromosome:
    def __init__(self, segments, name, length, origin_id=None):
        # Segments is a list of (start, end, origin_id)
        self.segments = segments
        self.name = name # '1', '2', ... 'X', 'Y'
        self.total_length = length
        # For pure chromosomes, if origin_id provided, ensure segments align
        if origin_id and not segments:
            self.segments = [(0, length, origin_id)]

    def length(self):
        return self.total_length

def create_pure_chromosome(name, origin):
    length = CHROMOSOME_LENGTHS[name]
    return Chromosome([(0, length, origin)], name, length)

def recombine(chr1, chr2):
    """
    Simulates meiosis with crossover between two homologous chromosomes.
    """
    # If distinct sex chromosomes (XY), no recombination in this simple model
    if chr1.name != chr2.name:
         # Standard XY non-recombining (ignoring PAR)
         return chr1 if random.random() < 0.5 else chr2

    length = chr1.length()
    # Number of crossovers depends on length
    expected_crossovers = length * CROSSOVER_RATE_FACTOR
    # Ensure at least one crossover often happens for large ones, but allow 0 for small
    num_crossovers = np.random.poisson(expected_crossovers)
    
    points = sorted([random.uniform(0, length) for _ in range(num_crossovers)])
    points = [0] + points + [length]

    new_segments = []
    
    # Pick starting chromosome randomly
    current_source = chr1 if random.random() < 0.5 else chr2
    other_source = chr2 if current_source == chr1 else chr1
    
    for i in range(len(points) - 1):
        start, end = points[i], points[i+1]
        
        # Extract slices from current_source
        chunk_start = start
        while chunk_start < end:
            # Find overlapping segment
            # Segment format: (seg_start, seg_end, origin)
            src_seg = next(s for s in current_source.segments if s[0] <= chunk_start < s[1])
            
            seg_end = min(end, src_seg[1])
            new_segments.append((chunk_start, seg_end, src_seg[2]))
            chunk_start = seg_end
            
        # Swap strands
        current_source, other_source = other_source, current_source

    # Merge adjacent
    merged = []
    if new_segments:
        curr_start, curr_end, curr_orig = new_segments[0]
        for next_start, next_end, next_orig in new_segments[1:]:
            # Floating point tolerance could be added, but exact match usually fine here
            if next_orig == curr_orig and abs(next_start - curr_end) < 1e-9:
                curr_end = next_end
            else:
                merged.append((curr_start, curr_end, curr_orig))
                curr_start, curr_end, curr_orig = next_start, next_end, next_orig
        merged.append((curr_start, curr_end, curr_orig))
    
    return Chromosome(merged, chr1.name, length)

class Person:
    def __init__(self, name, gender, genome=None):
        self.name = name
        self.gender = gender 
        self.genome = genome if genome else [] # List of pairs (chrA, chrB)

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
    
    # Map dad's gamete by chromosome name to find the right one (since X/Y order might vary or match index)
    # Actually, the lists are ordered 1..22, then Sex.
    
    for i in range(len(mom_gamete)):
        m_chr = mom_gamete[i]
        
        # Dad's corresponding chromosome
        # If it's the sex chromosome (last one)
        if i == 22: # Index 22 is the 23rd pair
            # Logic to force gender
            dad_sex_pair = father.genome[-1]
            
            # Identify which is X and which is Y in Dad's pair
            dad_X = dad_sex_pair[0] if dad_sex_pair[0].name == 'X' else dad_sex_pair[1]
            dad_Y = dad_sex_pair[0] if dad_sex_pair[0].name == 'Y' else dad_sex_pair[1]
            
            if child_gender == 'M':
                f_chr = dad_Y # Son gets Y
            else:
                f_chr = dad_X # Daughter gets X
                # Ideally, X can recombine if it's XX (Mother), but Dad is XY.
                # Dad's X passes unrecombined (pseudoautosomal regions ignored)
        else:
            f_chr = dad_gamete[i]
            
        child_genome.append([m_chr, f_chr])
        
    return Person(child_name, child_gender, child_genome)

# --- Setup Generation 1 (Grandparents) ---
def make_founder(name, gender, color):
    genome = []
    # Autosomes 1-22
    for i in range(1, 23):
        c1 = create_pure_chromosome(i, color)
        c2 = create_pure_chromosome(i, color)
        genome.append([c1, c2])
        
    # Sex Chromosomes
    if gender == 'M':
        cx = create_pure_chromosome('X', color)
        cy = create_pure_chromosome('Y', color)
        genome.append([cx, cy])
    else:
        c1 = create_pure_chromosome('X', color)
        c2 = create_pure_chromosome('X', color)
        genome.append([c1, c2])
        
    return Person(name, gender, genome)

MGM = make_founder("MGM", 'F', 'MGM')
MGF = make_founder("MGF", 'M', 'MGF')
PGM = make_founder("PGM", 'F', 'PGM')
PGF = make_founder("PGF", 'M', 'PGF')

# --- Setup Generation 2 (Parents) ---
# We use reproduce to generate them, but since grandparents are pure, 
# recombination just results in pure chromosomes of that color anyway.
# This correctly sets up the "Mother has one Red set, one Orange set" logic.
Mother = reproduce(MGM, MGF, 'F', "Mother")
Father = reproduce(PGM, PGF, 'M', "Father")

# --- Setup Generation 3 (Siblings) ---
Son = reproduce(Mother, Father, 'M', "Son")
Daughter = reproduce(Mother, Father, 'F', "Daughter")


# --- Plotting ---
# Wider figure to fit 23 pairs
fig, ax = plt.subplots(figsize=(20, 12))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

positions = {
    'MGM': (15, 88), 'MGF': (35, 88), 'PGM': (65, 88), 'PGF': (85, 88),
    'Mother': (25, 60), 'Father': (75, 60),
    'Daughter': (25, 25), 'Son': (75, 25) # Spread them out more so they can fit full width
}

def draw_person_genome(ax, person, center_x, center_y, width_limit=30):
    # width_limit is approx percent of canvas width this person can occupy
    
    ax.text(center_x, center_y + 8, person.name, ha='center', fontsize=14, fontweight='bold')
    
    num_pairs = 23
    # Spacing calculations
    # We have 23 pairs.
    # available width approx width_limit units.
    pair_spacing = width_limit / num_pairs 
    chr_width = pair_spacing * 0.4 # Bars are 40% of the slot
    
    scale_y = 0.05 # Scale height. Max len 250 -> 12.5 units tall
    
    start_x = center_x - (width_limit / 2) + (pair_spacing/2)
    
    for pair_idx, pair in enumerate(person.genome):
        pair_x = start_x + (pair_idx * pair_spacing)
        
        # Label the chromosome number below the pair
        chr_name = str(pair_idx + 1) if pair_idx < 22 else "X/Y"
        # Only label every few if crowded, or small font
        if pair_idx % 2 == 0 or pair_idx == 22:
             ax.text(pair_x, center_y - 8, chr_name, ha='center', fontsize=6, color='#555')
        
        for ch_idx, chrom in enumerate(pair):
            offset_x = - (chr_width * 0.6) if ch_idx == 0 else (chr_width * 0.6)
            x = pair_x + offset_x - (chr_width/2)
            
            # Align tops
            base_y = center_y 
            
            # Draw segments
            for seg in chrom.segments:
                seg_start, seg_end, origin = seg
                # Invert Y so 0 is top
                # chrom coordinate 0 is top.
                rect_y = base_y - (seg_end * scale_y)
                rect_h = (seg_end - seg_start) * scale_y
                rect_y = base_y - (seg_end * scale_y) # Bottom of rect
                
                # Matplotlib Rectangle is (x, y), w, h. y is bottom left.
                # If we want 0 at top:
                # Top of chromosome is at base_y.
                # Segment 0-10 starts at base_y, goes to base_y - 10*scale.
                # So rect y = base_y - seg_end*scale
                # rect h = (seg_end - seg_start)*scale
                
                color = GP_COLORS[origin]
                border = GP_BORDERS[origin]
                
                hatch = ''
                if origin in ['MGF', 'PGF']:
                    hatch = '////' if origin == 'MGF' else '||||'
                
                rect = patches.Rectangle((x, rect_y), chr_width, rect_h, 
                                         linewidth=0.3, edgecolor=border, facecolor=color, hatch=hatch)
                ax.add_patch(rect)
                
            # If it's Y, maybe label it specifically if pair is XY
            if chrom.name == 'Y':
                ax.text(x + chr_width/2, base_y - (chrom.length()*scale_y) - 2, "Y", ha='center', fontsize=5)

# Update positions to give more space
# Grandparents need to be squeezed or we make the figure huge.
# Let's give everyone equal width approx 40% of a quadrant.
# Actually, the siblings are the main event.
# Let's put Grandparents in top row, Parents middle, Children bottom.
# But 23 pairs is WIDE.
# Let's arrange them: 
# Row 1: MGM (0-25), MGF (25-50), PGM (50-75), PGF (75-100) -> 25% width each.
# Row 2: Mother (Left Half), Father (Right Half)
# Row 3: Daughter (Left Half), Son (Right Half)

# Clear axis and restart layout logic
ax.clear()
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Layout Configuration
row1_y = 85
row2_y = 55
row3_y = 20

# Width available for each person's genome block
gp_width = 22 # 4 people fit in 100
p_width = 40  # 2 people fit in 100
c_width = 45  # 2 people fit in 100

draw_person_genome(ax, MGM, 12.5, row1_y, gp_width)
draw_person_genome(ax, MGF, 37.5, row1_y, gp_width)
draw_person_genome(ax, PGM, 62.5, row1_y, gp_width)
draw_person_genome(ax, PGF, 87.5, row1_y, gp_width)

draw_person_genome(ax, Mother, 25, row2_y, p_width)
draw_person_genome(ax, Father, 75, row2_y, p_width)

draw_person_genome(ax, Daughter, 25, row3_y, c_width)
draw_person_genome(ax, Son, 75, row3_y, c_width)


# Connections (Symbolic lines)
def connect(x1, y1, x2, y2):
    ax.plot([x1, x2], [y1, y2], color='black', alpha=0.2, linewidth=1)

# GP to P
connect(12.5, row1_y-10, 25, row2_y+10) # MGM -> Mom
connect(37.5, row1_y-10, 25, row2_y+10) # MGF -> Mom
connect(62.5, row1_y-10, 75, row2_y+10) # PGM -> Dad
connect(87.5, row1_y-10, 75, row2_y+10) # PGF -> Dad

# P to C
connect(25, row2_y-10, 25, row3_y+10) # Mom -> Daughter
connect(75, row2_y-10, 25, row3_y+10) # Dad -> Daughter
connect(25, row2_y-10, 75, row3_y+10) # Mom -> Son
connect(75, row2_y-10, 75, row3_y+10) # Dad -> Son


# Legend
legend_elements = [
    patches.Patch(facecolor=GP_COLORS['MGM'], edgecolor=GP_BORDERS['MGM'], label='Maternal GM'),
    patches.Patch(facecolor=GP_COLORS['MGF'], edgecolor=GP_BORDERS['MGF'], hatch='////', label='Maternal GF'),
    patches.Patch(facecolor=GP_COLORS['PGM'], edgecolor=GP_BORDERS['PGM'], label='Paternal GM'),
    patches.Patch(facecolor=GP_COLORS['PGF'], edgecolor=GP_BORDERS['PGF'], hatch='||||', label='Paternal GF'),
]
ax.legend(handles=legend_elements, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0))

plt.suptitle("Human Chromosome Inheritance: 23 Pairs", fontsize=20, y=0.98)
plt.tight_layout()
plt.subplots_adjust(top=0.92) # Make room for title
plt.savefig("chromosome_inheritance_full.png", dpi=200)
print("Saved to chromosome_inheritance_full.png")
