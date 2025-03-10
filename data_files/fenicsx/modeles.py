from ufl import Identity, tr, grad, sym, inner, dx


class ModeleMecaniqueLineaireElastiqueIsotrope:
    def __init__(self, E: float, nu: float, hypothese_2D: str):
        # Stockage (et conversion) des paramètres élastiques
        self.mu = E / (2 * (1 + nu))
        self.la = E * nu / ((1 + nu) * (1 - 2 * nu))
        # Calcul des paramètres dépendants de l'hypothèse 2D
        # (contrainte de Kolosov k, coefficient formule d'Irwin Ep)
        match hypothese_2D:
            case "deformations_planes":
                self.k = 3 - 4 * nu
                self.Ep = E / (1 - nu**2)
            case "contraintes_planes":
                self.k = (3 - nu) / (1 + nu)
                self.Ep = E
                # Mise à jour coefficient de Lamé
                # NOTE: https://comet-fenics.readthedocs.io/en/latest/demo/elasticity/2D_elasticity.py.html
                self.la = 2 * self.mu * self.la / (self.la + 2 * self.mu)

    def eps(self, u):
        return sym(grad(u))

    def sig(self, u):
        mu, la = self.mu, self.la
        return 2 * mu * self.eps(u) + la * tr(self.eps(u)) * Identity(len(u))

    def energie_elastique(self, u):
        return 1 / 2 * inner(self.eps(u), self.sig(u)) * dx
