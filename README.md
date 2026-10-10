# Miroirs d'Archimède : simulation par lancer de rayons

Les soldats de Syracuse auraient-ils pu incendier la flotte romaine en renvoyant le soleil avec leurs boucliers de bronze ? Ce dépôt répond par le calcul : un moteur de lancer de rayons Monte-Carlo en C++ suit la lumière du soleil jusqu'aux miroirs puis jusqu'à la coque, et mesure le flux reçu en W/m², à comparer au seuil d'inflammation du bois.

![Carte de flux](figures/fig2_carte_flux.png)

## Ce que fait le programme

Pour chaque miroir, on tire des milliers de rayons. Chaque rayon :

1. part d'un point aléatoire du miroir ;
2. arrive du soleil avec une direction tirée dans le **disque solaire** (demi-angle 4,65 mrad), et non d'un point ;
3. porte une puissance `DNI × aire × cos(incidence) / nombre de rayons` ;
4. se réfléchit sur une surface légèrement ondulée (défauts de polissage, écart-type 2 mrad) ;
5. frappe le plan de la coque, où sa puissance multipliée par la réflectivité est déposée dans une grille de 5 cm.

Chaque miroir est orienté selon la bissectrice entre le soleil et le point visé, puis **dévié d'une erreur aléatoire** qui représente le soldat qui tient son bouclier à la main. Le pic de flux est lu sur une moyenne glissante de 15 cm × 15 cm pour ne pas confondre un point chaud réel avec du bruit statistique.

La position du soleil vient des formules d'astronomie de position (latitude de Syracuse, déclinaison d'été), l'éclairement direct du modèle de ciel clair de Meinel.

## Validation

Avant de croire un résultat, on vérifie le code sur un cas dont on connaît la réponse exacte. Pour un miroir plan parfait, le flux au centre de la tache vaut

  E = ρ · DNI · cos γ · F

où F est la fraction du disque solaire visible à travers le miroir depuis la coque. Tant que le miroir paraît plus grand que le soleil (D < D\* = a cos θ / 2α ≈ 81 m pour un miroir de 0,8 m), F = 1 et le flux est constant ; au-delà, il décroît en 1/D². Sur 17 distances de 5 à 300 m, la simulation colle à la théorie avec un écart moyen de 0,7 % (2,8 % au pire, ce qui correspond au bruit Monte-Carlo).

![Validation](figures/fig1_validation.png)

Cette figure contient déjà le résultat physique central : **un miroir plan ne concentre pas la lumière**. Il ne peut jamais envoyer plus de ρ × DNI, soit environ 0,55 kW/m² pour du bronze. Pour atteindre 20 kW/m², il faut donc superposer au minimum une quarantaine de taches, parfaitement alignées.

## Résultats

Seuil retenu : **20 kW/m²** pour l'auto-inflammation du bois avec flamme, et 10 kW/m² pour l'inflammation par incandescence (valeurs expérimentales de Boonmee, Université du Maryland).

**L'ennemi n'est pas le nombre de boucliers, c'est la précision de visée.**

![Miroirs requis](figures/fig3_miroirs_requis.png)

| Erreur de visée (par axe) | Navire à 20 m | à 50 m | à 100 m |
|---|---|---|---|
| 0° (visée parfaite) | ~39 | ~47 | ~109 |
| ±0,25° | ~47 | ~123 | ~470 |
| ±0,5° | ~86 | ~460 | plus de 500 |
| ±1° | ~380 | plus de 500 | plus de 500 |

Une erreur de 0,5° sur l'orientation d'un miroir dévie le faisceau de 1°, soit 87 cm à 50 m : les taches ne se superposent plus et le flux plafonne bien en dessous du seuil, quel que soit le nombre de soldats.

![Pic vs N](figures/fig4_pic_vs_N.png)

**L'heure compte aussi.** Le port de Syracuse s'ouvre à l'est. Le matin, soleil et navires sont du même côté : les miroirs travaillent presque face au soleil. L'après-midi, il faut renvoyer la lumière vers l'arrière, l'incidence devient rasante, la section du faisceau s'effondre. Vers 17 h solaires, même 500 miroirs ne suffisent plus.

![Heure](figures/fig5_heure.png)

## Hypothèses et limites

Ce modèle répond à une question précise (quel flux atteint la coque ?) et laisse de côté :

- l'**ombrage** et le **blocage** entre miroirs voisins (optimiste) ;
- la **forme réelle des boucliers**, convexes et ronds, ici supposés plans et carrés (optimiste : un miroir convexe disperse la lumière) ;
- le **mouvement du navire** et le temps nécessaire pour réorienter tous les miroirs ;
- le **temps d'exposition** : le seuil de 20 kW/m² doit être tenu plusieurs dizaines de secondes, il faudrait un modèle thermique de la planche ;
- l'**absorptivité** du bois dans le spectre solaire, différente de celle mesurée sous lampe infrarouge ;
- la réflectivité du bronze antique (0,6 ici), mal connue.

Presque toutes ces simplifications avantagent Archimède : les nombres de miroirs donnés sont des **minimums**.

## Compiler et lancer

```bash
cmake -S . -B build
cmake --build build
./build/archimede validate          # un miroir parfait contre la théorie (~30 s)
./build/archimede map N=130 D=50    # carte de flux pour un scénario
./build/archimede sweep             # distances x erreurs de visée (~2 min)
./build/archimede day D=50          # effet de l'heure (~1 min)
python3 scripts/plot_results.py     # figures dans figures/
```

Tout paramètre se règle en ligne de commande, par exemple `aim=0.5` (erreur de visée en degrés), `rho=0.8`, `cote=1.0`, `heure=15`, `pente=5`. La liste complète est dans `src/config.hpp`.

## Organisation du code

| Fichier | Rôle |
|---|---|
| `src/vec3.hpp` | vecteurs 3D, produit scalaire et vectoriel, loi de la réflexion |
| `src/config.hpp` | tous les paramètres et leur lecture en ligne de commande |
| `src/soleil.hpp` | position du soleil et éclairement direct |
| `src/simulation.hpp/.cpp` | miroirs, tirage des rayons, grille de flux |
| `src/main.cpp` | les quatre modes et l'écriture des CSV |
| `scripts/plot_results.py` | figures |

Les balayages utilisent OpenMP s'il est disponible : chaque simulation est indépendante, elles sont réparties sur tous les coeurs du processeur.

## Pistes pour aller plus loin

- Ajouter l'ombrage entre miroirs (tester si le rayon incident ou réfléchi croise un autre miroir).
- Modéliser des boucliers ronds et légèrement convexes.
- Coupler la carte de flux à un modèle de conduction 1D dans la planche pour obtenir un temps d'inflammation.
- Comparer aux essais expérimentaux réels (MIT, 2005 ; MythBusters).

## Références

- N. Boonmee, *Theoretical and Experimental Study of Autoignition of Wood*, thèse, University of Maryland, 2004. https://drum.lib.umd.edu/handle/1903/1691
- A. B. Meinel et M. P. Meinel, *Applied Solar Energy*, Addison-Wesley, 1976 (modèle de ciel clair).
