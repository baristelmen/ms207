from ufl import Identity, tr, grad, sym, inner, dx


class ModeleMecaniqueLineaireElastiqueIsotrope:
    def __init__(self, E: float, nu: float):
        # Stockage (et conversion) des paramètres élastiques
        self.la = E * nu / ((1 + nu) * (1 - nu))
        self.mu = E / (2 * (1 + nu))
        # Calcul de la constante de Kosolov (déformations planes)
        self.k = 3 - 4 * nu

    def eps(self, u):
        return 1 / 2 * sym(grad(u))

    def sig(self, u):
        mu, la = self.mu, self.la
        return 2 * mu * self.eps(u) + la * tr(self.eps(u)) * Identity(len(u))

    def energie_elastique(self, u):
        return 1 / 2 * inner(self.eps(u), self.sig(u)) * dx
