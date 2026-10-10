// vec3.hpp
// Petit vecteur 3D pour le lancer de rayons.
// Repère utilisé partout : x vers l'est (la mer), y vers le nord, z vers le haut.
#pragma once

#include <cmath>

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec3 operator-(const Vec3& b) const { return {x - b.x, y - b.y, z - b.z}; }
    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3 operator*(double k) const { return {k * x, k * y, k * z}; }
};

inline Vec3 operator*(double k, const Vec3& a) { return a * k; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double norme(const Vec3& a) { return std::sqrt(dot(a, a)); }

inline Vec3 normaliser(const Vec3& a) { return a * (1.0 / norme(a)); }

// Réflexion spéculaire d'une direction d sur une surface de normale n (unitaire).
// Loi de Snell-Descartes pour la réflexion : r = d - 2 (d.n) n
inline Vec3 reflechir(const Vec3& d, const Vec3& n) { return d - 2.0 * dot(d, n) * n; }

// Construit deux vecteurs u, v tels que (u, v, n) soit une base orthonormée directe.
// u est choisi horizontal quand c'est possible (bord "horizontal" du miroir).
inline void baseOrthonormee(const Vec3& n, Vec3& u, Vec3& v) {
    Vec3 haut(0.0, 0.0, 1.0);
    Vec3 c = cross(n, haut);
    if (norme(c) < 1e-9) c = cross(n, Vec3(0.0, 1.0, 0.0));  // n vertical : cas dégénéré
    u = normaliser(c);
    v = cross(n, u);
}
