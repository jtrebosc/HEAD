#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <sstream>
#include <stdexcept>

#include "simpson_rw.hpp"

// Fonction pour lire un fichier .spe ou .fid
bool readSimpsonFile(const std::string& filename, SSpectrum& spec) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << std::endl;
        return false;
    }

    std::string line;

    // Lire la première ligne (doit contenir "SIMP")
    if (!std::getline(file, line) || line != "SIMP") {
        std::cerr << "Erreur : Le fichier ne commence pas par 'SIMP'." << std::endl;
        file.close();
        return false;
    }

    // Lire les paramètres
    while (std::getline(file, line) && line != "DATA") {
        std::istringstream iss(line);
        std::string param;
        if (std::getline(iss, param, '=')) {
            if (param == "NP") {
                iss >> spec.NP;
            } else if (param == "SW") {
                iss >> spec.SW;
            } else if (param == "REF") { // REF is the position of 0 ppm in Hz with respect to center of spectrum (carrier)
                iss >> spec.REF;
            } else if (param == "SFREQ") {
                iss >> spec.SFREQ;
            } else if (param == "NI") {
                iss >> spec.NI;
            } else if (param == "SW1") {
                iss >> spec.SW1;
            } else if (param == "REF1") {
                iss >> spec.REF1;
            } else if (param == "SFREQ1") {
                iss >> spec.SFREQ1;
            } else if (param == "TYPE") {
                iss >> spec.TYPE;
            }
        }
    }

    // Lire les données complexes
    spec.ComplexData.resize(spec.NI, std::vector<std::complex<double>>(spec.NP));
    for (int i = 0; i < spec.NI; ++i) {
        for (int j = 0; j < spec.NP; ++j) {
            double real, imag;
            if (!(file >> real >> imag)) {
                std::cerr << "Erreur : Impossible de lire les données complexes." << std::endl;
                file.close();
                return false;
            }
            spec.ComplexData[i][j] = std::complex<double>(real, imag);
        }
    }

    // Vérifier la fin du fichier
    std::string endMarker;
    if (!(file >> endMarker) || endMarker != "END") {
        std::cerr << "Erreur : Le fichier ne se termine pas par 'END'." << std::endl;
        file.close();
        return false;
    }

    file.close();
    return true;
}

// Fonction pour écrire un fichier .spe ou .fid
bool writeSimpsonFile(const std::string& filename, const SSpectrum& spec) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible de créer le fichier " << filename << std::endl;
        return false;
    }

    // Écrire l'en-tête
    file << "SIMP\n";
    file << "NP=" << spec.NP << "\n";
    file << "SW=" << spec.SW << "\n";
    file << "REF=" << spec.REF << "\n";
    file << "SFREQ=" << spec.SFREQ << "\n";
    file << "NI=" << spec.NI << "\n";
    file << "SW1=" << spec.SW1 << "\n";
    file << "REF1=" << spec.REF1 << "\n";
    file << "SFREQ1=" << spec.SFREQ1 << "\n";
    file << "TYPE=" << spec.TYPE << "\n";
    file << "DATA\n";

    // Écrire les données complexes
    for (const auto& row : spec.ComplexData) {
        for (const auto& value : row) {
            file << value.real() << " " << value.imag() << "\n";
        }
    }

    // Écrire la fin du fichier
    file << "END\n";
    file.close();
    return true;
}
/*

// Exemple d'utilisation
int main() {
    SSpectrum spec;

    // Exemple de lecture
    if (readSimpsonFile("examples/original_spectrum.spe", spec)) {
        std::cout << "Fichier lu avec succès !" << std::endl;
        std::cout << "NP: " << spec.NP << ", NI: " << spec.NI << std::endl;
        std::cout << "SW: " << spec.SW << ", REF: " << spec.REF << std::endl;
        std::cout << "SFREQ: " << spec.SFREQ << std::endl;
        std::cout << "SW1: " << spec.SW1 << ", REF1: " << spec.REF1 << std::endl;
        std::cout << "SFREQ1: " << spec.SFREQ1 << std::endl;

        // Afficher les premières données complexes
        std::cout << "Premières données complexes : " << std::endl;
        for (size_t i = 0; i < std::min(spec.ComplexData.size(), static_cast<size_t>(2)); ++i) {
            for (size_t j = 0; j < std::min(spec.ComplexData[i].size(), static_cast<size_t>(2)); ++j) {
                std::cout << spec.ComplexData[i][j] << " ";
            }
            std::cout << std::endl;
        }
    } else {
        std::cerr << "Échec de la lecture du fichier." << std::endl;
    }

    // Exemple d'écriture
    SSpectrum newspec; //= {4, 2, 10000.0, 100.0, 400.13, 20000.0, 200.0, 101.03, "SPE"};
    newspec.NP = 4;
    newspec.NI = 2;
    newspec.SW = 10000.;
    newspec.SW1 = 20000.;
    newspec.TYPE = "SPE";
    newspec.ComplexData = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}},
        {{1.1, 1.0}, {2.1, 1.0}, {3.1, 1.0}, {4.1, 1.0}}
    };

    if (writeSimpsonFile("new_spectrum.spe", newspec)) {
        std::cout << "Fichier écrit avec succès !" << std::endl;
    } else {
        std::cerr << "Échec de l'écriture du fichier." << std::endl;
    }

    return 0;
}
*/
