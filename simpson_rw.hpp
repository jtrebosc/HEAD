#include <complex>
#include <vector>
#include <string>

// Structure pour stocker les paramètres du spectre
struct SSpectrum { // Simpson spectrum or FID structure
    int NP=1;         // Nombre de points dans la dimension directe
    int NI=1;         // Nombre de points dans la dimension indirecte
    double SW=1;      // Fenêtre spectrale (dimension directe)
    double REF=0;     // Valeur en ppm du centre du spectre (dimension directe)
    double SFREQ=1;   // Fréquence absolue du centre du spectre (dimension directe)
    double SW1=1;     // Fenêtre spectrale (dimension indirecte)
    double REF1=0;    // Valeur en ppm du centre du spectre (dimension indirecte)
    double SFREQ1=1;  // Fréquence absolue du centre du spectre (dimension indirecte)
    std::string TYPE="SPE";  // Type du spectre: spectre fréquentiel (SPE) ou Free Induction Decay temporel (FID)
    std::vector<std::vector<std::complex<double>>> ComplexData;
};

// Structure pour stocker les données complexes

// Fonction pour lire un fichier .spe ou .fid
bool readSimpsonFile(const std::string& filename, SSpectrum& spec) ;
