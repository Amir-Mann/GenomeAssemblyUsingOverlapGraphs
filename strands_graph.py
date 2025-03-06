"""
Module for constructing and using an overlap graph for genome assembly.

This module includes classes and functions for building a trie for efficient suffix-prefix matching
and constructing a graph of genome read overlaps. It also provides utilities for visualizing the trie and
timing specific operations.
"""

import math
import time
from contextlib import contextmanager

timings = {}


@contextmanager
def timer(key: str):
    """
    Context manager to measure execution time of a code block.
    """
    start = time.time()
    yield
    if key not in timings:
        timings[key] = 0
    timings[key] += time.time() - start


class Vertex:
    """
    Represents a vertex in the overlap graph, holding a genome read sequence and connections to other vertices.
    """
    def __init__(self, sequence: str) -> None:
        """
        Initializes a Vertex with the given sequence.
        """
        self.sequence: str = sequence
        self.connected_vertices: list['Vertex'] = []


class TrieNode:
    """
    Represents a node in a trie used for fast suffix-prefix matching of genome reads.
    """
    def __init__(self, extra: Vertex = None, count: int = 0) -> None:
        """
        Initializes a TrieNode.
        """
        self.A: 'TrieNode' = None
        self.C: 'TrieNode' = None
        self.G: 'TrieNode' = None
        self.T: 'TrieNode' = None
        self.count: int = count
        self.extra: Vertex | None = extra

    def str_aux(self):
        """
        Helper recursive function to generate a string representation of the trie for visualization.
        """
        sub_tree = []
        lines = 0
        if self.A is not None:
            A_sub_tree, A_lines = self.A.str_aux()
            sub_tree.append("-\033[30;105mA\033[0m" + A_sub_tree[0])
            sub_tree.extend([" \\" + line for line in A_sub_tree[1:]])
            lines += A_lines
        if self.C is not None:
            C_sub_tree, C_lines = self.C.str_aux()
            sub_tree.append("-\033[30;105mC\033[0m" + C_sub_tree[0])
            sub_tree.extend([" \\" + line for line in C_sub_tree[1:]])
            lines += C_lines
        if self.G is not None:
            G_sub_tree, G_lines = self.G.str_aux()
            sub_tree.append("-\033[30;105mG\033[0m" + G_sub_tree[0])
            sub_tree.extend([" \\" + line for line in G_sub_tree[1:]])
            lines += G_lines
        if self.T is not None:
            T_sub_tree, T_lines = self.T.str_aux()
            sub_tree.append("-\033[30;105mT\033[0m" + T_sub_tree[0])
            sub_tree.extend([" \\" + line for line in T_sub_tree[1:]])
            lines += T_lines
        if lines == 0:
            lines += 1
            sub_tree = [""]
        return sub_tree, lines
    
    def __str__(self) -> str:
        """
        Returns a string representation of the trie.
        """
        return "\n".join(self.str_aux()[0])

    def insert_strand(self, strand: str, extra: Vertex | None = None, new_extra: bool = False) -> Vertex | None:
        """
        Inserts a strand into the trie and associates it with a Vertex.
        """
        self.count += 1
        if self.extra is None:
            self.extra = extra if extra is not None else Vertex(strand)
            new_extra = True
        
        if not strand:
            return self.extra if new_extra else None
        
        first_char = strand[0]
        next_node: TrieNode | None = getattr(self, first_char)
        if next_node is None:  # Opening recursion for faster runtimes.
            current = self
            extra = self.extra if new_extra else Vertex(strand)
            for c in strand:
                next_node = TrieNode(extra, 1)
                setattr(current, c, next_node)
                current = next_node
            return extra
        
        return next_node.insert_strand(strand[1:], self.extra if new_extra else None, new_extra=new_extra)

    def search(self, strand: str, allow_mis_matches: int = 0) -> Vertex | None:
        """
        Searches for a strand in the trie, allowing for a specified number of mismatches.
        """
        node = self
        for i, char in enumerate(strand):
            next_node: TrieNode | None = getattr(node, char)
            if next_node is None:
                if not allow_mis_matches:
                    return None
                for char in "ACGT":
                    extra = node.search(char + strand[i+1:], allow_mis_matches=allow_mis_matches - 1)
                    if extra is not None:
                        return extra
                return None
            node = next_node
        return node.extra


class Graph:
    """
    Represents an overlap graph where nodes are genome read vertices and edges represent overlaps.
    """
    def __init__(self) -> None:
        """
        Initializes an empty Graph.
        """
        self.vertices: dict[Vertex, None] = {}

    def load_from_strands(self, strands: list[str], allow_mis_matches: int) -> None:
        """
        Constructs the overlap graph from a list of DNA strands using a trie for efficient matching.
        """
        trie = TrieNode()
        
        for strand in strands:
            with timer("insert_strand"):
                vertex = trie.insert_strand(strand)
            if vertex is not None:
                vertex.sequence = strand
                self.vertices[vertex] = None
        
        for vertex in self.vertices.keys():
            for i in range(1, len(vertex.sequence) - 2 * int(math.log(len(strands), 4))):
                suffix_strand = vertex.sequence[i:]
                with timer("search"):
                    match = trie.search(suffix_strand, allow_mis_matches)
                if match and match != vertex:
                    vertex.connected_vertices.append((match, len(vertex.sequence) - i))
                    break

    def get_sequenced_result(self, min_overlap: int = 7, in_coming_links_min: int = 40) -> str:
        """
        Traverses the overlap graph to generate the assembled genome sequence.
        """
        incoming_links_values = {v: 0 for v in self.vertices.keys()}
        
        for vertex in self.vertices.keys():
            for edge_end, match_length in vertex.connected_vertices:
                incoming_links_values[edge_end] += match_length

        max_seq = ""

        for origin, sum_incoming_links in incoming_links_values.items():
            if sum_incoming_links > in_coming_links_min:
                continue
                
            current = origin
            sequence = current.sequence
            i = 0
            while current.connected_vertices and i < len(self.vertices):
                next_vertex, overlap_size = current.connected_vertices[0]
                if overlap_size < min_overlap:
                    break
                
                sequence += next_vertex.sequence[overlap_size:]
                current = next_vertex
                i += 1
            
            if len(sequence) > len(max_seq):
                max_seq = sequence

        return max_seq
