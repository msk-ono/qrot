#ifndef QROT_PARSER_H
#define QROT_PARSER_H

#include <memory>
#include <string>
#include <vector>

#include "qrot/boost.h"

namespace qrot {
enum class TokenKind { RESERVED, NUM, END };
struct Token {
    TokenKind kind;
    std::string_view str;
};
enum class NodeKind { ADD, SUB, MUL, DIV, NUM };
struct Node {
    NodeKind kind;
    const Token* token;
    const Node* l;
    const Node* r;
};
/**
 * @brief Parses strings representing arithmetic expressions involving the basic operations
 * (+ - * /) and supports the use of 'pi' as a keyword.
 * @details This class is implemented with reference to https://www.sigbus.info/compilerbook.
 *
 * Definition of expressions.
 * expr    = mul ("+" mul | "-" mul)*
 * mul     = unary ("*" unary | "/" unary)*
 * unary   = ("+" | "-")? primary
 * primary = num | "(" expr ")"
 */
class AST {
public:
    static AST Parse(const std::string& s);

    AST() = default;
    AST(const AST&) = delete;  // Delete copy
    AST(AST&&) = default;
    AST& operator=(const AST&) = delete;  // Delete copy
    AST& operator=(AST&&) = default;

    const Node* Root() const { return root_; }
    Float Value() const;

private:
    const Node* root_;
    std::unique_ptr<std::string> str_;
    std::unique_ptr<std::vector<Token>> tokens_;
    std::unique_ptr<std::vector<Node>> nodes_;
};
}  // namespace qrot

#endif  // QROT_PARSER_H
