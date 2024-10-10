import itertools
import networkx as nx
import matplotlib.pyplot as plt


# Función para calcular el solapamiento entre dos secuencias
def calcular_solapamiento(seq1, seq2):
    max_overlap = 0
    for i in range(1, min(len(seq1), len(seq2)) + 1):
        if seq1[-i:] == seq2[:i]:
            max_overlap = i
    return max_overlap


# Función para crear un grafo de solapamientos entre secuencias
def crear_grafo_solapamiento(secuencias):
    G = nx.DiGraph()
    n = len(secuencias)

    for i in range(n):
        for j in range(n):
            if i != j:
                overlap = calcular_solapamiento(secuencias[i], secuencias[j])
                if overlap > 0:
                    G.add_edge(
                        i, j, weight=-overlap
                    )  # Peso negativo para minimizar superstring

    return G


# Función para encontrar el camino hamiltoniano óptimo
def encontrar_camino_hamiltoniano(G, secuencias):
    n = len(secuencias)
    mejor_superstring = None
    mejor_longitud = float("inf")

    # Generar todas las posibles permutaciones (caminos hamiltonianos)
    for permutacion in itertools.permutations(range(n)):
        superstring = secuencias[permutacion[0]]

        for i in range(1, n):
            u = permutacion[i - 1]
            v = permutacion[i]
            overlap = calcular_solapamiento(secuencias[u], secuencias[v])
            superstring += secuencias[v][overlap:]

        if len(superstring) < mejor_longitud:
            mejor_superstring = superstring
            mejor_longitud = len(superstring)

    return mejor_superstring


# Función para visualizar el grafo de solapamientos
def visualizar_grafo(G, secuencias):
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G, "weight")

    nx.draw(
        G,
        pos,
        with_labels=True,
        labels={i: secuencias[i] for i in range(len(secuencias))},
        node_color="lightblue",
        node_size=2000,
        font_size=10,
        font_weight="bold",
        arrows=True,
    )
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title("Grafo de Solapamientos entre Secuencias")
    plt.show()


# Secuencias de prueba
secuencias = [
    "ATCCGTTGAAGCCGCGGGC",
    "TTAACTCGAGG",
    "TTAAGTACTGCCCG",
    "ATCTGTGTCGGG",
    "CGACTCCCGACACA",
    "CACAGATCCGTTGAAGCCGCGGG",
    "CTCGAGTTAAGTA",
    "CGCGGGCAGTACTT",
]

# Crear el grafo de solapamientos
G = crear_grafo_solapamiento(secuencias)

# Visualizar el grafo
visualizar_grafo(G, secuencias)

# Encontrar el camino hamiltoniano óptimo
superstring_optimo = encontrar_camino_hamiltoniano(G, secuencias)
print("Secuencia consenso mínima (superstring):", superstring_optimo)
