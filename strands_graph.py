

class Vertex:
    def __init__(self, sequence: str) -> None:
        self.sequence: str = sequence
        self.connected_vertices: list['Vertex'] = []


class TrieNode:
    def __init__(self) -> None:
        self.A: 'TrieNode' = None
        self.C: 'TrieNode' = None
        self.G: 'TrieNode' = None
        self.T: 'TrieNode' = None
        self.count: int = 0
        self.extra: Vertex | None = None

    def insert_strand(self, strand: str, extra: Vertex | None = None) -> Vertex | None:
        self.count += 1
        new_extra = False
        if self.extra is None:
            self.extra = extra if extra is not None else Vertex(strand)
            new_extra = True
        
        if not strand:
            return self.extra if new_extra else None
        
        first_char = strand[0]
        next_node: TrieNode | None = getattr(self, first_char)
        if next_node is None:
            next_node = TrieNode()
            setattr(self, first_char, next_node)
        
        return next_node.insert_strand(strand[1:], self.extra)

    def search(self, strand: str) -> Vertex | None:
        node = self
        for char in strand:
            next_node: TrieNode | None = getattr(node, char)
            if next_node is None:
                return node.extra
            node = next_node
        return node.extra


class Graph:
    def __init__(self) -> None:
        self.vertices: dict[Vertex, None] = {}

    def load_from_strands(self, strands: list[str]) -> None:
        trie = TrieNode()
        
        for strand in strands:
            vertex = trie.insert_strand(strand)
            if vertex is not None:
                vertex.sequence = strand
                self.vertices[vertex] = None
        
        for vertex in self.vertices.keys():
            reverse_strand = vertex.sequence[::-1]
            best_match = trie.search(reverse_strand)
            if best_match and best_match != vertex:
                vertex.connected_vertices.append(best_match)

    def get_sequenced_result(self) -> str:
        incoming_links = {v: False for v in self.vertices.keys()}
        
        for vertex in self.vertices.keys():
            for connected in vertex.connected_vertices:
                incoming_links[connected] = True
        
        for origin, really_an_origin in incoming_links.items():
            if not really_an_origin:
                break
        else:
            print("Error no origin")
            return ""
        
        current = origin
        sequence = current.sequence
        
        while current.connected_vertices:
            next_vertex = current.connected_vertices[0]
            
            max_overlap = 0
            for i in range(1, min(len(current.sequence), len(next_vertex.sequence)) + 1):
                if current.sequence[-i:] == next_vertex.sequence[:i]:
                    max_overlap = i
                else:
                    break
            
            sequence += next_vertex.sequence[max_overlap:]
            current = next_vertex
        
        return sequence
