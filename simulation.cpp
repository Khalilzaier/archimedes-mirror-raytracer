// simulation.cpp
#include "simulation.hpp"

#include <algorithm>
#include <cmath>

// =====================================================================
//  GrilleFlux
// =====================================================================

GrilleFlux::GrilleFlux(double demiLargeur, double pas)
    : demi_(demiLargeur), pas_(pas), n_(static_cast<int>(std::round(2.0 * demiLargeur / pas))),
      brut_(n_ * n_, 0.0), lisse_(n_ * n_, 0.0) {}

bool GrilleFlux::deposer(double y, double z, double puissance) {
    int i = static_cast<int>(std::floor((y + demi_) / pas_));
    int j = static_cast<int>(std::floor((z + demi_) / pas_));
    if (i < 0 || j < 0 || i >= n_ || j >= n_) return false;

    brut_[index(i, j)] += puissance;
    total_ += puissance;

    // La maille (i, j) contribue à la moyenne 3x3 de ses 9 voisines.
    // Comme les valeurs ne font qu'augmenter, le maximum se met à jour au fil de l'eau.
    for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
            int a = i + di, b = j + dj;
            if (a < 0 || b < 0 || a >= n_ || b >= n_) continue;
            double& cell = lisse_[index(a, b)];
            cell += puissance;
            if (cell > picLisse_) picLisse_ = cell;
        }
    }
    return true;
}

double GrilleFlux::flux(int i, int j) const { return brut_[index(i, j)] / (pas_ * pas_); }

double GrilleFlux::fluxLisse(int i, int j) const { return lisse_[index(i, j)] / (9.0 * pas_ * pas_); }

double GrilleFlux::surfaceAuDessus(double seuil) const {
    int compte = 0;
    for (int j = 0; j < n_; ++j)
        for (int i = 0; i < n_; ++i)
            if (fluxLisse(i, j) >= seuil) ++compte;
    return compte * pas_ * pas_;
}

// =====================================================================
//  Simulation
// =====================================================================

Simulation::Simulation(const Config& c, const Vec3& soleil, double dni, unsigned graine)
    : cfg_(c), soleil_(soleil), dni_(dni), grille_(c.demi_largeur, c.pas_grille), rng_(graine) {}

Vec3 Simulation::positionMiroir(int k) const {
    int col = k % cfg_.colonnes;
    int rangee = k / cfg_.colonnes;
    double y = (col - 0.5 * (cfg_.colonnes - 1)) * cfg_.pas_lateral;
    double x = -rangee * cfg_.recul_rangee;
    double z = cfg_.hauteur_quai + rangee * cfg_.montee_rangee;
    return {x, y, z};
}

// Tire une direction uniformément répartie (en angle solide) dans le disque solaire.
// Le soleil n'est pas un point : c'est un disque de 0,53 degré. C'est précisément ce
// qui empêche un miroir plan de concentrer la lumière en un point.
Vec3 Simulation::tirerDirectionSoleil(const Vec3& u, const Vec3& v) {
    double cosMax = std::cos(cfg_.alpha_soleil);
    double cosT = 1.0 - uniforme_(rng_) * (1.0 - cosMax);
    double sinT = std::sqrt(1.0 - cosT * cosT);
    double phi = 2.0 * 3.14159265358979323846 * uniforme_(rng_);
    return cosT * soleil_ + sinT * (std::cos(phi) * u + std::sin(phi) * v);
}

// Incline légèrement une normale : petits angles gaussiens d'écart-type sigma sur deux axes.
Vec3 Simulation::perturber(const Vec3& n, const Vec3& u, const Vec3& v, double sigma) {
    if (sigma <= 0.0) return n;
    return normaliser(n + sigma * gauss_(rng_) * u + sigma * gauss_(rng_) * v);
}

void Simulation::ajouterMiroir() {
    Miroir m;
    m.centre = positionMiroir(nombreMiroirs());

    // Orientation idéale : la normale est la bissectrice entre la direction du soleil
    // et la direction de la cible (c'est la loi de la réflexion à l'envers).
    Vec3 versCible = normaliser(pointVise() - m.centre);
    Vec3 nIdeale = normaliser(soleil_ + versCible);

    // Le soldat oriente son bouclier à la main : erreur aléatoire, fixée pour ce miroir.
    Vec3 u0, v0;
    baseOrthonormee(nIdeale, u0, v0);
    m.normale = perturber(nIdeale, u0, v0, cfg_.sigma_visee_rad());
    baseOrthonormee(m.normale, m.u, m.v);
    miroirs_.push_back(m);

    // Repère du disque solaire, pour tirer les directions des rayons incidents.
    Vec3 su, sv;
    baseOrthonormee(soleil_, su, sv);

    // Plan de la cible : vertical, de normale tournée vers le rivage.
    const Vec3 P0 = pointVise();
    const Vec3 nCible(-1.0, 0.0, 0.0);
    const Vec3 axeY(0.0, 1.0, 0.0), axeZ(0.0, 0.0, 1.0);

    const double a = cfg_.cote_miroir;
    const double aire = a * a;
    const int R = cfg_.rayons;
    const double sigmaPente = cfg_.sigma_pente_rad();

    double puissanceMiroir = 0.0;

    for (int r = 0; r < R; ++r) {
        // 1. Point d'impact aléatoire sur la surface du miroir
        Vec3 p = m.centre + (uniforme_(rng_) - 0.5) * a * m.u + (uniforme_(rng_) - 0.5) * a * m.v;

        // 2. Direction d'arrivée tirée dans le disque solaire
        Vec3 versSoleil = tirerDirectionSoleil(su, sv);
        double cosI = dot(m.normale, versSoleil);
        if (cosI <= 0.0) continue;  // le soleil éclaire le dos du miroir

        // 3. Puissance portée par ce rayon : DNI x surface x cos(incidence) / nombre de rayons
        double w = dni_ * aire * cosI / R;
        puissanceMiroir += w;

        // 4. Réflexion sur une surface légèrement ondulée (défauts de polissage)
        Vec3 nLocale = perturber(m.normale, m.u, m.v, sigmaPente);
        Vec3 dir = reflechir(-versSoleil, nLocale);

        // 5. Intersection avec le plan de la coque
        double denom = dot(dir, nCible);
        if (std::fabs(denom) < 1e-12) continue;
        double t = dot(P0 - p, nCible) / denom;
        if (t <= 0.0) continue;  // le rayon part dans l'autre sens
        Vec3 impact = p + t * dir;

        // 6. Dépôt de l'énergie réfléchie (fraction rho) dans la grille
        Vec3 rel = impact - P0;
        grille_.deposer(dot(rel, axeY), dot(rel, axeZ), cfg_.reflectivite * w);
    }

    puissanceMiroirs_ += puissanceMiroir;
    sommeIncidence_ += std::acos(std::clamp(dot(m.normale, soleil_), -1.0, 1.0));
}

double Simulation::angleIncidenceMoyen() const {
    return miroirs_.empty() ? 0.0 : sommeIncidence_ / miroirs_.size();
}
