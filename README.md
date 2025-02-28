# Projet de Programmation Parallèle

_« Architecture matérielle et logicielle des ordinateurs »_

## Instructions

### Objectif

Écrire un programme « multi-threadé » et « vectoriel » en C, puis comprendre les divers aspects liés à son fonctionnement correct ainsi qu’à ses performances.

### Remarque

Le programme doit être écrit, **en individuel**, en langage C, en considérant **Pthread et AVX**.

### Exercice 1

Écrire la fonction

```c
double dist(float *U, float *V, int n);
```

qui calcule

$$
dist(u,v,n)=\sum_{i=0}^{n-1}\sqrt(\frac{u_i^2 + v_i^2}{1+(u_iv_i)^2})
$$

pour deux vecteurs U et V de n réels en simple précision.

### Exercice 2

Écrire une version vectorielle (AVX) de la fonction précédente en supposant que les adresses de U et V sont alignées et que n est un multiple de 8 (**ind** : utiliser les **intrinsics** pour la racine carrée **vus en TPs**).

```c
double vect_dist(float *U, float *V, int n);
```

### Exercice 3

Écrire une version plus générale de la fonction `vect_dist` qui fonctionne **sans contrainte d’alignement** (**ind** : _loadu et storeu_) mais en supposant toujours que n est un multiple de 8.

```c
double vect_dist_gen(float *U, float *V, int n);
```

### Exercice 4

Proposer une version **multithreadée** de la fonction `dist()`

```c
void distPar(float *U, float *V, int n, int nb_threads, int mode);
```

qui considère une exécution **multi-threadée** avec `nb_threads` threads. On peut par exemple supposer que chaque thread traitera un **bloc** des vecteurs U et V (**partitionnement statique**).

Le paramètre `mode` indique le type de calcul (**0 : scalaire | 1 : vectoriel**).

Penser à créer **une structure pour échanger les données du thread**.

**NB :**

- Utilisez **le mutex pour gérer la somme** et **appelez `vect_dist()` ou `vect_dist_gen()`**.
- Si `p` est un pointeur sur un vecteur dont on s’intéresse aux valeurs `p[k]` avec k allant de i à j−1, on peut utiliser le pointeur `q = p + i` et considérer `q[k]` avec k allant de 0 et (j−i)−1.

### Exercice 5

Écrire la fonction principale du projet.

Le `main()` doit :

a) Créer et **initialiser les vecteurs** U et V avec des **valeurs aléatoires dans `]0,1]`**.

b) **Calculer et afficher** (en prenant `n=1024²`) :

- `dist(U, V, n)`,
- **Temps d’exécution séquentielle**,
- **Temps d’exécution parallèle**,
- **Accélération**.

## Rendu

Voici quelques remarques sur le rendu proposé:

- Separation du code dans plusieurs fichiers pour facilité la lisibilité (compilation via le makefile proposé). Les fichiers d'interet principal sont `operation.c` et `main.c` pour les implémentations des différentes fonctions proposées.
- Le calcul vectoriel étant réalisé sur des floats et non des doubles, mon implementation souffre d'erreurs d'approximation... Je fourni donc également la fonction `vect_dist_double` avec les caractéristiques suivantes et son benchmark.

```c
double vect_dist_double(double *U, double *V, int n);
```
