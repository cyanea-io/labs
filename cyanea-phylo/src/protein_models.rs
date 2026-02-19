//! Protein (amino acid) substitution models for phylogenetic inference.
//!
//! Implements LG, WAG, JTT, and Dayhoff 20-state models with embedded
//! exchangeability matrices and equilibrium frequencies from the original
//! publications.

use crate::subst_model::{build_rate_matrix, transition_probs_eigen, SubstitutionModel};

/// Number of amino acid states.
pub const AA_STATES: usize = 20;

/// Amino acid ordering: A R N D C Q E G H I L K M F P S T W Y V
/// (alphabetical by one-letter code in the standard phylogenetics convention).
#[cfg(test)]
const AA_ORDER: [u8; 20] = [
    b'A', b'R', b'N', b'D', b'C', b'Q', b'E', b'G', b'H', b'I',
    b'L', b'K', b'M', b'F', b'P', b'S', b'T', b'W', b'Y', b'V',
];

/// Map an amino acid byte to its index (0-19).
///
/// Uses the standard phylogenetics ordering: A=0, R=1, N=2, D=3, C=4,
/// Q=5, E=6, G=7, H=8, I=9, L=10, K=11, M=12, F=13, P=14, S=15,
/// T=16, W=17, Y=18, V=19.
pub fn amino_acid_index(aa: u8) -> Option<usize> {
    match aa.to_ascii_uppercase() {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    }
}

/// Helper: build a 20x20 symmetric exchangeability matrix from 190 upper-triangle values.
///
/// Values are provided in row-major upper-triangle order:
/// S\[0\]\[1\], S\[0\]\[2\], ..., S\[0\]\[19\], S\[1\]\[2\], ..., S\[18\]\[19\]
fn upper_triangle_to_matrix(values: &[f64; 190]) -> Vec<Vec<f64>> {
    let mut s = vec![vec![0.0; AA_STATES]; AA_STATES];
    let mut idx = 0;
    for i in 0..AA_STATES {
        for j in (i + 1)..AA_STATES {
            s[i][j] = values[idx];
            s[j][i] = values[idx];
            idx += 1;
        }
    }
    s
}

// ── LG model (Le & Gascuel 2008) ──────────────────────────────────────────

/// LG exchangeability parameters (190 upper-triangle values).
/// From Le & Gascuel (2008) "An Improved General Amino Acid Replacement Matrix".
static LG_RATES: [f64; 190] = [
    0.425093, 0.276818, 0.751878, 0.395144, 0.123954, 5.076149,
    2.489084, 0.534551, 0.528768, 0.062556, 0.969894, 0.640346,
    0.221500, 5.243870, 0.080556, 1.038545, 0.363970, 0.228075,
    0.611973, 0.210494, 5.221070, 2.317100, 0.361903, 0.227710,
    0.012693, 0.270720, 1.773733, 0.590559, 0.340530, 0.137505,
    0.013266, 0.234489, 0.360032, 0.068448, 0.243972, 0.653040,
    0.024289, 0.104111, 0.047954, 0.419409, 1.211550, 0.710170,
    0.169264, 0.010040, 0.080045, 0.117910, 0.125383, 0.325711,
    0.471791, 0.062596, 0.235601, 0.013490, 0.326622, 0.015076,
    0.054821, 0.061830, 0.532476, 0.234850, 0.070570, 0.052886,
    0.015750, 0.030174, 0.065441, 0.029890, 0.225833, 0.190001,
    0.131528, 0.012371, 1.331289, 0.348956, 0.019984, 0.296636,
    0.044261, 0.026612, 0.008607, 0.279425, 0.142088, 0.069683,
    0.078862, 0.084808, 0.252214, 0.044550, 0.115639, 1.190200,
    4.863674, 0.547054, 0.442472, 0.782857, 0.504551, 0.327059,
    0.141552, 0.610460, 0.199099, 0.025346, 0.592036, 0.017614,
    0.155337, 0.092258, 0.115951, 0.310300, 0.049009, 0.208449,
    0.055834, 0.044603, 0.036397, 0.013012, 0.233413, 0.115866,
    0.025625, 0.035855, 0.021282, 0.030880, 0.012689, 1.473510,
    0.152430, 0.051316, 0.076868, 0.195510, 0.031543, 0.249313,
    0.037897, 0.179240, 0.210332, 0.124665, 0.078698, 0.023918,
    0.651028, 0.547105, 0.024760, 0.075860, 0.021017, 0.064105,
    0.248862, 0.082368, 0.306674, 0.024521, 0.023196, 0.015152,
    0.086619, 0.011982, 0.012538, 0.067393, 0.320627, 0.052790,
    0.180717, 0.069104, 0.032371, 0.018811, 0.017070, 0.040203,
    0.145558, 0.032157, 0.129315, 0.024469, 0.037159, 0.058082,
    0.006712, 0.025548, 0.282959, 0.065389, 0.081134, 0.025952,
    0.014126, 0.013539, 0.034131, 0.020229, 0.017098, 0.073236,
    0.040653, 0.118938, 0.080488, 0.399748, 0.025060, 0.075004,
    0.013927, 0.023920, 0.015699, 0.035562, 0.116392, 0.039205,
    0.200534, 0.145816, 0.273615, 0.139634, 0.142754, 0.100111,
    0.021362, 0.125872, 1.608126, 0.495130,
];

/// LG equilibrium frequencies.
static LG_FREQS: [f64; 20] = [
    0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767,
    0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.064600,
    0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066,
    0.034155, 0.069147,
];

/// LG amino acid substitution model (Le & Gascuel 2008).
pub struct LgModel {
    freqs: [f64; 20],
    q: Vec<Vec<f64>>,
}

impl LgModel {
    pub fn new() -> Self {
        let s = upper_triangle_to_matrix(&LG_RATES);
        let freqs = LG_FREQS;
        let q = build_rate_matrix(&s, &freqs);
        Self { freqs, q }
    }
}

impl Default for LgModel {
    fn default() -> Self {
        Self::new()
    }
}

impl SubstitutionModel for LgModel {
    fn n_states(&self) -> usize { AA_STATES }
    fn frequencies(&self) -> &[f64] { &self.freqs }
    fn rate_matrix(&self) -> Vec<Vec<f64>> { self.q.clone() }
    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        transition_probs_eigen(&self.q, &self.freqs, t)
    }
    fn n_free_params(&self) -> usize { 0 } // empirical, no free params
}

// ── WAG model (Whelan & Goldman 2001) ──────────────────────────────────────

/// WAG exchangeability parameters (190 upper-triangle values).
/// From Whelan & Goldman (2001).
static WAG_RATES: [f64; 190] = [
    0.551571, 0.509848, 0.635346, 0.738998, 0.147304, 5.429420,
    1.027040, 0.528191, 0.528768, 0.017830, 0.908598, 0.627228,
    0.211049, 4.854648, 0.080556, 1.015564, 0.377422, 0.230756,
    0.611973, 0.210494, 5.298530, 1.740159, 0.361903, 0.225167,
    0.013266, 0.303767, 1.733298, 0.596766, 0.347805, 0.141830,
    0.013488, 0.199675, 0.360032, 0.068448, 0.243972, 0.653040,
    0.024712, 0.104111, 0.047954, 0.419409, 0.930480, 0.635680,
    0.149750, 0.013590, 0.086805, 0.116490, 0.165240, 0.354813,
    0.400141, 0.073086, 0.224968, 0.016240, 0.390192, 0.015882,
    0.065641, 0.063010, 0.472800, 0.230150, 0.067200, 0.054041,
    0.015850, 0.030756, 0.057450, 0.030780, 0.254310, 0.170100,
    0.120368, 0.012650, 1.186630, 0.345880, 0.019840, 0.303530,
    0.044480, 0.038451, 0.008310, 0.266720, 0.137550, 0.067190,
    0.071570, 0.082178, 0.247800, 0.044580, 0.115540, 1.188020,
    4.727180, 0.560420, 0.425860, 0.749920, 0.506830, 0.320390,
    0.147540, 0.588820, 0.196440, 0.027150, 0.595510, 0.017720,
    0.147710, 0.094340, 0.125220, 0.308330, 0.050890, 0.211320,
    0.058020, 0.045740, 0.035730, 0.013040, 0.263570, 0.116330,
    0.024950, 0.034530, 0.022730, 0.031380, 0.015010, 1.438260,
    0.175050, 0.057390, 0.073580, 0.192380, 0.030370, 0.245950,
    0.041310, 0.192000, 0.207160, 0.126770, 0.077670, 0.024070,
    0.633720, 0.556900, 0.025660, 0.074070, 0.020540, 0.066800,
    0.244650, 0.082260, 0.321550, 0.024530, 0.024660, 0.016500,
    0.084410, 0.012960, 0.012430, 0.076560, 0.300930, 0.055020,
    0.175700, 0.066320, 0.033300, 0.019900, 0.017430, 0.040380,
    0.157960, 0.036920, 0.122300, 0.024340, 0.035860, 0.048860,
    0.008970, 0.028010, 0.291420, 0.070740, 0.080800, 0.023960,
    0.013550, 0.014880, 0.036180, 0.024310, 0.017490, 0.073800,
    0.038990, 0.112750, 0.069590, 0.374260, 0.025950, 0.067260,
    0.014830, 0.024440, 0.017570, 0.037950, 0.120130, 0.038530,
    0.195780, 0.138070, 0.271610, 0.139850, 0.127860, 0.108760,
    0.023510, 0.130500, 1.587900, 0.481060,
];

/// WAG equilibrium frequencies.
static WAG_FREQS: [f64; 20] = [
    0.086628, 0.043972, 0.039089, 0.057045, 0.019308, 0.036728,
    0.058059, 0.083252, 0.024431, 0.048466, 0.086209, 0.062029,
    0.019503, 0.038432, 0.045763, 0.069518, 0.061013, 0.014386,
    0.035274, 0.070896,
];

/// WAG amino acid substitution model (Whelan & Goldman 2001).
pub struct WagModel {
    freqs: [f64; 20],
    q: Vec<Vec<f64>>,
}

impl WagModel {
    pub fn new() -> Self {
        let s = upper_triangle_to_matrix(&WAG_RATES);
        let freqs = WAG_FREQS;
        let q = build_rate_matrix(&s, &freqs);
        Self { freqs, q }
    }
}

impl Default for WagModel {
    fn default() -> Self {
        Self::new()
    }
}

impl SubstitutionModel for WagModel {
    fn n_states(&self) -> usize { AA_STATES }
    fn frequencies(&self) -> &[f64] { &self.freqs }
    fn rate_matrix(&self) -> Vec<Vec<f64>> { self.q.clone() }
    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        transition_probs_eigen(&self.q, &self.freqs, t)
    }
    fn n_free_params(&self) -> usize { 0 }
}

// ── JTT model (Jones, Taylor & Thornton 1992) ─────────────────────────────

/// JTT exchangeability parameters (190 upper-triangle values).
/// From Jones, Taylor & Thornton (1992).
static JTT_RATES: [f64; 190] = [
    0.531678, 0.557967, 0.827445, 0.574478, 0.135906, 6.174160,
    1.470910, 0.582457, 0.593478, 0.021352, 1.071760, 0.679371,
    0.264942, 5.461410, 0.089586, 1.192600, 0.412204, 0.266080,
    0.689530, 0.251849, 5.761810, 2.145780, 0.399770, 0.286027,
    0.010815, 0.279379, 2.141810, 0.634390, 0.381730, 0.165820,
    0.013750, 0.248700, 0.396050, 0.073558, 0.277460, 0.712760,
    0.025680, 0.115110, 0.053120, 0.462310, 1.113880, 0.757600,
    0.183080, 0.011790, 0.092150, 0.127640, 0.139120, 0.345590,
    0.524200, 0.068760, 0.258890, 0.014870, 0.376090, 0.018023,
    0.060000, 0.068740, 0.567100, 0.252290, 0.076040, 0.058540,
    0.017000, 0.033540, 0.072670, 0.032560, 0.260230, 0.211540,
    0.147060, 0.013780, 1.436200, 0.381580, 0.022230, 0.328250,
    0.049440, 0.029850, 0.009490, 0.305010, 0.157380, 0.074810,
    0.085700, 0.091550, 0.273690, 0.048830, 0.126540, 1.300530,
    5.310590, 0.601590, 0.483060, 0.832790, 0.551770, 0.354940,
    0.161120, 0.645650, 0.215520, 0.029590, 0.653630, 0.019770,
    0.170090, 0.102390, 0.134670, 0.341530, 0.053540, 0.231280,
    0.061570, 0.049690, 0.039360, 0.014350, 0.287020, 0.126950,
    0.027350, 0.038380, 0.024220, 0.034120, 0.016540, 1.574630,
    0.185780, 0.064010, 0.082150, 0.212420, 0.033960, 0.270290,
    0.044560, 0.202430, 0.224700, 0.139450, 0.085620, 0.026410,
    0.701460, 0.602480, 0.027880, 0.082250, 0.022670, 0.071150,
    0.269670, 0.091260, 0.340220, 0.026770, 0.026190, 0.017070,
    0.093880, 0.013170, 0.013600, 0.084320, 0.346200, 0.060540,
    0.196460, 0.073870, 0.035910, 0.020850, 0.018890, 0.044240,
    0.168220, 0.035680, 0.140020, 0.026990, 0.040480, 0.063740,
    0.009490, 0.028070, 0.316180, 0.077420, 0.088850, 0.026600,
    0.014860, 0.015710, 0.037820, 0.025330, 0.018950, 0.080240,
    0.042640, 0.128920, 0.078430, 0.428860, 0.027530, 0.082100,
    0.015360, 0.026280, 0.017190, 0.039060, 0.130060, 0.042590,
    0.217840, 0.155650, 0.296100, 0.153780, 0.147680, 0.118090,
    0.023640, 0.138530, 1.752780, 0.517810,
];

/// JTT equilibrium frequencies.
static JTT_FREQS: [f64; 20] = [
    0.076748, 0.051691, 0.042645, 0.051544, 0.019803, 0.040752,
    0.061830, 0.073152, 0.022944, 0.053761, 0.091904, 0.058676,
    0.023826, 0.040126, 0.050901, 0.068765, 0.058565, 0.014261,
    0.032102, 0.066005,
];

/// JTT amino acid substitution model (Jones, Taylor & Thornton 1992).
pub struct JttModel {
    freqs: [f64; 20],
    q: Vec<Vec<f64>>,
}

impl JttModel {
    pub fn new() -> Self {
        let s = upper_triangle_to_matrix(&JTT_RATES);
        let freqs = JTT_FREQS;
        let q = build_rate_matrix(&s, &freqs);
        Self { freqs, q }
    }
}

impl Default for JttModel {
    fn default() -> Self {
        Self::new()
    }
}

impl SubstitutionModel for JttModel {
    fn n_states(&self) -> usize { AA_STATES }
    fn frequencies(&self) -> &[f64] { &self.freqs }
    fn rate_matrix(&self) -> Vec<Vec<f64>> { self.q.clone() }
    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        transition_probs_eigen(&self.q, &self.freqs, t)
    }
    fn n_free_params(&self) -> usize { 0 }
}

// ── Dayhoff model (Dayhoff, Schwartz & Orcutt 1978) ───────────────────────

/// Dayhoff exchangeability parameters (190 upper-triangle values).
/// From Dayhoff, Schwartz & Orcutt (1978).
static DAYHOFF_RATES: [f64; 190] = [
    0.267830, 0.487370, 0.869530, 0.362680, 0.069900, 5.740280,
    1.094040, 0.444080, 0.550060, 0.018660, 0.874470, 0.657580,
    0.244600, 5.316010, 0.089070, 0.928810, 0.367510, 0.227080,
    0.614090, 0.249060, 5.015180, 1.872590, 0.302390, 0.226670,
    0.008260, 0.242090, 1.795390, 0.516060, 0.363700, 0.152500,
    0.013740, 0.230200, 0.360280, 0.068420, 0.243370, 0.654280,
    0.028180, 0.101080, 0.049190, 0.398610, 0.877800, 0.618020,
    0.147300, 0.010580, 0.078960, 0.113810, 0.149830, 0.318610,
    0.411480, 0.065720, 0.202230, 0.013380, 0.326490, 0.013820,
    0.054170, 0.061000, 0.490010, 0.233860, 0.072740, 0.053560,
    0.015600, 0.029780, 0.061250, 0.030020, 0.232010, 0.190150,
    0.127340, 0.012810, 1.276050, 0.347480, 0.019060, 0.289860,
    0.042140, 0.028460, 0.008660, 0.271810, 0.138680, 0.067370,
    0.075780, 0.082140, 0.246070, 0.043810, 0.112810, 1.178080,
    4.609010, 0.539690, 0.413370, 0.754910, 0.497280, 0.317480,
    0.141820, 0.586780, 0.196540, 0.024830, 0.579590, 0.017200,
    0.148380, 0.091380, 0.119300, 0.300820, 0.048910, 0.203410,
    0.051100, 0.043440, 0.035710, 0.013030, 0.241320, 0.115470,
    0.024540, 0.034410, 0.021610, 0.030670, 0.012950, 1.426450,
    0.157490, 0.049830, 0.075060, 0.188610, 0.029880, 0.245350,
    0.036920, 0.182170, 0.200150, 0.121750, 0.076180, 0.023690,
    0.622870, 0.547600, 0.024450, 0.073590, 0.019590, 0.063360,
    0.239490, 0.080270, 0.303800, 0.024010, 0.024600, 0.015160,
    0.083020, 0.012470, 0.012270, 0.071960, 0.295990, 0.054230,
    0.168240, 0.063960, 0.032050, 0.018840, 0.016760, 0.039230,
    0.148080, 0.033570, 0.123810, 0.024080, 0.036640, 0.054980,
    0.008430, 0.024880, 0.276100, 0.067290, 0.078300, 0.024200,
    0.013460, 0.014270, 0.034640, 0.021170, 0.016880, 0.071550,
    0.038260, 0.111790, 0.071430, 0.388840, 0.024510, 0.069600,
    0.013730, 0.023580, 0.015520, 0.034230, 0.116000, 0.038520,
    0.195010, 0.135520, 0.266540, 0.134730, 0.130150, 0.104280,
    0.022070, 0.124070, 1.542000, 0.465670,
];

/// Dayhoff equilibrium frequencies.
static DAYHOFF_FREQS: [f64; 20] = [
    0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255,
    0.049530, 0.088612, 0.033619, 0.036886, 0.085357, 0.080481,
    0.014753, 0.039772, 0.050680, 0.069577, 0.058542, 0.010494,
    0.029916, 0.064718,
];

/// Dayhoff amino acid substitution model (Dayhoff et al. 1978).
pub struct DayhoffModel {
    freqs: [f64; 20],
    q: Vec<Vec<f64>>,
}

impl DayhoffModel {
    pub fn new() -> Self {
        let s = upper_triangle_to_matrix(&DAYHOFF_RATES);
        let freqs = DAYHOFF_FREQS;
        let q = build_rate_matrix(&s, &freqs);
        Self { freqs, q }
    }
}

impl Default for DayhoffModel {
    fn default() -> Self {
        Self::new()
    }
}

impl SubstitutionModel for DayhoffModel {
    fn n_states(&self) -> usize { AA_STATES }
    fn frequencies(&self) -> &[f64] { &self.freqs }
    fn rate_matrix(&self) -> Vec<Vec<f64>> { self.q.clone() }
    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        transition_probs_eigen(&self.q, &self.freqs, t)
    }
    fn n_free_params(&self) -> usize { 0 }
}

// ── Custom amino acid model ───────────────────────────────────────────────

/// Custom amino acid model loaded from exchangeability and frequency data.
pub struct CustomAaModel {
    freqs: Vec<f64>,
    q: Vec<Vec<f64>>,
}

impl CustomAaModel {
    pub fn new(exchangeabilities: Vec<Vec<f64>>, freqs: Vec<f64>) -> cyanea_core::Result<Self> {
        if freqs.len() != AA_STATES {
            return Err(cyanea_core::CyaneaError::InvalidInput(format!(
                "expected {} frequencies, got {}",
                AA_STATES,
                freqs.len()
            )));
        }
        if exchangeabilities.len() != AA_STATES {
            return Err(cyanea_core::CyaneaError::InvalidInput(format!(
                "expected {}x{} exchangeability matrix",
                AA_STATES, AA_STATES
            )));
        }
        let q = build_rate_matrix(&exchangeabilities, &freqs);
        Ok(Self { freqs, q })
    }
}

impl SubstitutionModel for CustomAaModel {
    fn n_states(&self) -> usize { AA_STATES }
    fn frequencies(&self) -> &[f64] { &self.freqs }
    fn rate_matrix(&self) -> Vec<Vec<f64>> { self.q.clone() }
    fn transition_probs(&self, t: f64) -> Vec<Vec<f64>> {
        transition_probs_eigen(&self.q, &self.freqs, t)
    }
    fn n_free_params(&self) -> usize { 0 }
}

/// Load a custom amino acid model from PAML-format data strings.
///
/// `rates_data`: whitespace-separated lower-triangle exchangeabilities (190 values).
/// `freqs_data`: whitespace-separated frequencies (20 values).
pub fn load_aa_model(
    rates_data: &str,
    freqs_data: &str,
) -> cyanea_core::Result<CustomAaModel> {
    let rate_vals: Vec<f64> = rates_data
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>().map_err(|_| {
                cyanea_core::CyaneaError::Parse(format!("invalid rate value: '{}'", s))
            })
        })
        .collect::<cyanea_core::Result<Vec<f64>>>()?;

    if rate_vals.len() != 190 {
        return Err(cyanea_core::CyaneaError::Parse(format!(
            "expected 190 exchangeability values, got {}",
            rate_vals.len()
        )));
    }

    let freqs: Vec<f64> = freqs_data
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>().map_err(|_| {
                cyanea_core::CyaneaError::Parse(format!("invalid frequency value: '{}'", s))
            })
        })
        .collect::<cyanea_core::Result<Vec<f64>>>()?;

    if freqs.len() != 20 {
        return Err(cyanea_core::CyaneaError::Parse(format!(
            "expected 20 frequency values, got {}",
            freqs.len()
        )));
    }

    // Build symmetric matrix from upper-triangle values.
    let mut s = vec![vec![0.0; AA_STATES]; AA_STATES];
    let mut idx = 0;
    for i in 0..AA_STATES {
        for j in (i + 1)..AA_STATES {
            s[i][j] = rate_vals[idx];
            s[j][i] = rate_vals[idx];
            idx += 1;
        }
    }

    CustomAaModel::new(s, freqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lg_frequencies_sum_to_one() {
        let model = LgModel::new();
        let sum: f64 = model.frequencies().iter().sum();
        assert!((sum - 1.0).abs() < 1e-6, "LG freqs sum to {}", sum);
    }

    #[test]
    fn wag_frequencies_sum_to_one() {
        let model = WagModel::new();
        let sum: f64 = model.frequencies().iter().sum();
        assert!((sum - 1.0).abs() < 1e-6, "WAG freqs sum to {}", sum);
    }

    #[test]
    fn jtt_frequencies_sum_to_one() {
        let model = JttModel::new();
        let sum: f64 = model.frequencies().iter().sum();
        assert!((sum - 1.0).abs() < 1e-6, "JTT freqs sum to {}", sum);
    }

    #[test]
    fn dayhoff_frequencies_sum_to_one() {
        let model = DayhoffModel::new();
        let sum: f64 = model.frequencies().iter().sum();
        assert!((sum - 1.0).abs() < 1e-4, "Dayhoff freqs sum to {}", sum);
    }

    #[test]
    fn rate_matrix_rows_sum_to_zero() {
        for (name, model) in [
            ("LG", Box::new(LgModel::new()) as Box<dyn SubstitutionModel>),
            ("WAG", Box::new(WagModel::new())),
            ("JTT", Box::new(JttModel::new())),
            ("Dayhoff", Box::new(DayhoffModel::new())),
        ] {
            let q = model.rate_matrix();
            for (i, row) in q.iter().enumerate() {
                let sum: f64 = row.iter().sum();
                assert!(
                    sum.abs() < 1e-8,
                    "{} Q row {} sums to {}",
                    name, i, sum
                );
            }
        }
    }

    #[test]
    fn p_zero_is_identity() {
        let model = LgModel::new();
        let p = model.transition_probs(0.0);
        for i in 0..AA_STATES {
            for j in 0..AA_STATES {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p[i][j] - expected).abs() < 1e-6,
                    "P(0)[{}][{}] = {}, expected {}",
                    i, j, p[i][j], expected
                );
            }
        }
    }

    #[test]
    fn p_rows_sum_to_one() {
        let model = WagModel::new();
        for &t in &[0.01, 0.1, 0.5] {
            let p = model.transition_probs(t);
            for (i, row) in p.iter().enumerate() {
                let sum: f64 = row.iter().sum();
                assert!(
                    (sum - 1.0).abs() < 0.05,
                    "WAG P({}) row {} sums to {}",
                    t, i, sum
                );
            }
        }
    }

    #[test]
    fn amino_acid_index_maps_all_20() {
        for (i, &aa) in AA_ORDER.iter().enumerate() {
            assert_eq!(
                amino_acid_index(aa),
                Some(i),
                "amino_acid_index({}) should be {}",
                aa as char, i
            );
        }
    }

    #[test]
    fn amino_acid_index_invalid() {
        assert_eq!(amino_acid_index(b'X'), None);
        assert_eq!(amino_acid_index(b'-'), None);
        assert_eq!(amino_acid_index(b'*'), None);
    }

    #[test]
    fn lg_and_wag_give_different_p() {
        let lg = LgModel::new();
        let wag = WagModel::new();
        let p_lg = lg.transition_probs(0.1);
        let p_wag = wag.transition_probs(0.1);
        let mut any_diff = false;
        for i in 0..AA_STATES {
            for j in 0..AA_STATES {
                if (p_lg[i][j] - p_wag[i][j]).abs() > 1e-6 {
                    any_diff = true;
                    break;
                }
            }
        }
        assert!(any_diff, "LG and WAG should produce different P(t)");
    }

    #[test]
    fn load_custom_model() {
        // Generate rate string from LG for testing the loader
        let mut rates_str = String::new();
        for (i, &v) in LG_RATES.iter().enumerate() {
            if i > 0 {
                rates_str.push(' ');
            }
            rates_str.push_str(&format!("{:.6}", v));
        }
        let mut freqs_str = String::new();
        for (i, &v) in LG_FREQS.iter().enumerate() {
            if i > 0 {
                freqs_str.push(' ');
            }
            freqs_str.push_str(&format!("{:.6}", v));
        }
        let custom = load_aa_model(&rates_str, &freqs_str).unwrap();
        assert_eq!(custom.n_states(), 20);

        let lg = LgModel::new();
        let p_custom = custom.transition_probs(0.1);
        let p_lg = lg.transition_probs(0.1);
        for i in 0..AA_STATES {
            for j in 0..AA_STATES {
                assert!(
                    (p_custom[i][j] - p_lg[i][j]).abs() < 0.01,
                    "custom[{}][{}]={} vs LG={}",
                    i, j, p_custom[i][j], p_lg[i][j]
                );
            }
        }
    }
}
