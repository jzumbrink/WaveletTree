#include <memory>
#include <iostream>

struct Node {
    int value;
    std::unique_ptr<Node> left;
    std::unique_ptr<Node> right;
};

// Rekursive Funktion
void print_tree(const Node& node, int depth = 0) {
    std::cout << std::string(depth, '-') << node.value << "\n";
    if (node.left) print_tree(*node.left, depth + 1);
    if (node.right) print_tree(*node.right, depth + 1);
}

int main() {
    auto root = std::make_unique<Node>();
    root->value = 10;
    root->left = std::make_unique<Node>();
    root->left->value = 5;
    root->right = std::make_unique<Node>();
    root->right->value = 15;

    print_tree(*root);  // ✅ Übergabe per Referenz
}
