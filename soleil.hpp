// soleil.hpp
// Position du soleil (modèle simple d'éphémérides) et éclairement direct.
#pragma once

#include <algorithm>
#include <cmath>

#include "config.hpp"
#include "vec3.hpp"

constexpr double PI = 3.14159265358979323846;
inline double rad(double deg) { return deg * PI / 180.0; }
inline double deg(double r) { return r * 180.0 / PI; }

// Vecteur unitaire pointant vers le soleil, dans le repère (est, nord, haut).
// H = angle horaire (15 degrés par heure, négatif le matin), delta = déclinaison,
// phi = latitude. Formules classiques de l'astronomie de position.
inline Vec3 directionSoleil(const Config& c) {
    if (c.soleil_el_deg > -90.0 && c.soleil_az_deg > -900.0) {
        double el = rad(c.soleil_el_deg), az = rad(c.soleil_az_deg);
        return {std::cos(el) * std::sin(az), std::cos(el) * std::cos(az), std::sin(el)};
    }
    double H = rad(15.0 * (c.heure_solaire - 12.0));
    double d = rad(c.declinaison_deg);
    double phi = rad(c.latitude_deg);
    Vec3 s(-std::cos(d) * std::sin(H),
           std::sin(d) * std::cos(phi) - std::cos(d) * std::cos(H) * std::sin(phi),
           std::sin(d) * std::sin(phi) + std::cos(d) * std::cos(H) * std::cos(phi));
    return normaliser(s);
}

// Éclairement direct normal (W/m2) par ciel clair, modèle empirique de Meinel :
// DNI = 1353 * 0.7^(AM^0.678), avec la masse d'air AM = 1 / sin(élévation).
inline double eclairementDirect(const Config& c, const Vec3& s) {
    if (c.dni > 0.0) return c.dni;
    if (s.z <= 0.0) return 0.0;  // soleil sous l'horizon
    double am = 1.0 / std::max(s.z, 0.02);
    return 1353.0 * std::pow(0.7, std::pow(am, 0.678));
}
