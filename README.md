# Pricing d'Options - Projet C++

Implémentation des trois principales méthodes de pricing pour options financières : Black-Scholes, arbre binomial et Monte Carlo.

## Compilation

```bash
make
./pricing
```

Ou directement :
```bash
g++ -std=c++20 -O3 main.cpp -o pricing -lm
```

## Ce qui a été implémenté

### Trois méthodes de pricing

**Black-Scholes** : La formule analytique classique. Rapide (< 1ms) et exacte pour les options européennes. J'ai aussi codé le calcul des grecques (delta, gamma, vega, etc).

**Arbre binomial** : Méthode de Cox-Ross-Rubinstein. Plus lent mais permet de pricer les options américaines avec exercice anticipé. Converge bien vers Black-Scholes avec suffisamment de steps.

**Monte Carlo** : Simulations stochastiques avec le générateur Mersenne Twister. Indispensable pour les options path-dependent comme les asiatiques. Un peu lent (~1s pour 100k simulations) mais très flexible.

### Types d'options

- Options vanille (call/put, européennes/américaines)
- Options asiatiques (payoff sur prix moyen)

### Architecture du code

Trois fichiers header principaux :

- `Options.h` : Classes pour les différents types d'options
- `Market.h` : Paramètres de marché et structures de résultats  
- `Pricers.h` : Les trois algorithmes de pricing

Le code utilise du polymorphisme classique (classes abstraites, héritage) et des smart pointers pour éviter les fuites mémoire.

## Exemples de résultats

Pour un call européen standard (S=100, K=100, T=1 an, r=5%, σ=20%) :

| Méthode | Prix | Temps |
|---------|------|-------|
| Black-Scholes | 10.4506 | 0.03 ms |
| Binomial (1000 steps) | 10.4486 | 16 ms |
| Monte Carlo (100k) | 10.4474 | 770 ms |

Les trois méthodes convergent bien (erreur < 0.5%).

## Structure du code

```
Options.h       Classes d'options (Option, VanillaOption, AsianOption)
Market.h        Données de marché et résultats
Pricers.h       Algorithmes (BlackScholesPricer, BinomialTreePricer, MonteCarloPricer)
main.cpp        Programme avec exemples interactifs
Makefile        Compilation
```

## Points techniques

J'ai fait attention à quelques trucs :

- Utilisation de `std::shared_ptr` et `std::unique_ptr` plutôt que des pointeurs bruts
- Mersenne Twister (`std::mt19937`) pour la génération de nombres aléatoires
- Optimisation de la compilation avec `-O3` et `-march=native`
- Exception handling pour les paramètres invalides
- Const-correctness partout

Pour les options américaines, l'arbre binomial vérifie à chaque nœud si l'exercice anticipé est optimal (max entre exercice et continuation).

Les options asiatiques nécessitent de stocker tout le path de prix, c'est pour ça qu'elles ne marchent qu'avec Monte Carlo.

## Quelques détails mathématiques

### Black-Scholes

```
C = S·e^(-qT)·N(d1) - K·e^(-rT)·N(d2)

d1 = [ln(S/K) + (r-q+σ²/2)T] / (σ√T)
d2 = d1 - σ√T
```

### Binomial CRR

```
u = e^(σ√Δt)
d = 1/u
p = (e^(r-q)Δt - d) / (u - d)
```

### Monte Carlo

Mouvement brownien géométrique :
```
S(t+Δt) = S(t) · exp[(r-q-σ²/2)Δt + σ√Δt·Z]
```
où Z suit une loi normale N(0,1).

## Améliorations possibles

- Paralléliser Monte Carlo avec OpenMP (facile à faire)
- Ajouter d'autres options exotiques (barrières, lookback...)
- Implémenter Longstaff-Schwartz pour les américaines avec MC
- Calculer la volatilité implicite (Newton-Raphson)

## Dépendances

Aucune bibliothèque externe. Standard library C++20 uniquement.

---

Projet réalisé par : Nada Mousteau, Maxime Hebert, Stanislas de Chezelles
