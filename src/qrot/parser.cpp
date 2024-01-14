#include "parser.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace qrot {
namespace {
bool IsDigits(char c) { return '0' <= c && c <= '9'; }
void Tokenize(const std::string& str, std::vector<Token>& container) {
    const auto sv = std::string_view(str);
    for (auto i = std::size_t{0}; i < str.size(); ++i) {
        if (str[i] == ' ' || str[i] == '\n') {
            // Do nothing
        } else if (str[i] == '+' || str[i] == '-' || str[i] == '*' || str[i] == '/' ||
                   str[i] == '(' || str[i] == ')') {
            container.emplace_back(TokenKind::RESERVED, sv.substr(i, 1));
        } else if (str[i] == 'p') {
            if (i + 1 < str.size() && str[i + 1] == 'i') {
                container.emplace_back(TokenKind::NUM, sv.substr(i++, 2));
            } else {
                throw std::runtime_error("Found unknown token: " + str.substr(i));
            }
        } else if (IsDigits(str[i])) {
            auto num_dots = 0;
            auto j = i;
            while (j < str.size() && (IsDigits(str[j]) || str[j] == '.')) {
                if (str[j] == '.') {
                    if (num_dots != 0) { throw std::runtime_error("Too many dots in digits"); }
                    num_dots++;
                }
                j++;
            }
            container.emplace_back(TokenKind::NUM, sv.substr(i, j - i));
            i = j - 1;
        } else {
            throw std::runtime_error(std::string("Found unknown token: ") + str[i]);
        }
    }
    container.emplace_back(TokenKind::END);
}
void GoToNextToken(Token** tokens) { (*tokens)++; }
const Node* ConstructExpr(Token** tokens, std::vector<Node>& container);
const Node* ConstructMul(Token** tokens, std::vector<Node>& container);
const Node* ConstructPrimary(Token** tokens, std::vector<Node>& container);
const Node* ConstructUnary(Token** tokens, std::vector<Node>& container);
const Node* ConstructExpr(Token** tokens, std::vector<Node>& container) {
    const auto* node = ConstructMul(tokens, container);
    while (true) {
        const auto* const token = *tokens;
        if (token->kind == TokenKind::END) { break; }
        const auto c = token->str[0];
        if (c == '+') {
            GoToNextToken(tokens);
            container.emplace_back(NodeKind::ADD, token, node, ConstructMul(tokens, container));
            node = &container.back();
        } else if (c == '-') {
            GoToNextToken(tokens);
            container.emplace_back(NodeKind::SUB, token, node, ConstructMul(tokens, container));
            node = &container.back();
        } else {
            // END or ')
            break;
        }
    }
    return node;
}
const Node* ConstructMul(Token** tokens, std::vector<Node>& container) {
    const auto* node = ConstructUnary(tokens, container);
    while (true) {
        const auto* const token = *tokens;
        if (token->kind == TokenKind::END) { break; }
        const auto c = token->str[0];
        if (c == '*') {
            GoToNextToken(tokens);
            container.emplace_back(NodeKind::MUL, token, node, ConstructUnary(tokens, container));
            node = &container.back();
        } else if (c == '/') {
            GoToNextToken(tokens);
            container.emplace_back(NodeKind::DIV, token, node, ConstructUnary(tokens, container));
            node = &container.back();
        } else {
            break;
        }
    }
    return node;
}
const Node* ConstructUnary(Token** tokens, std::vector<Node>& container) {
    const auto* const token = *tokens;
    if (token->kind == TokenKind::END) {
        throw std::runtime_error(
            "Expected unary expression but the actual is the end of the string");
    }
    const auto c = token->str[0];
    if (c == '+') {
        GoToNextToken(tokens);
        return ConstructPrimary(tokens, container);
    } else if (c == '-') {
        GoToNextToken(tokens);
        container.emplace_back(NodeKind::SUB, token, ConstructPrimary(tokens, container), nullptr);
        return &container.back();
    }
    return ConstructPrimary(tokens, container);
}
const Node* ConstructPrimary(Token** tokens, std::vector<Node>& container) {
    const auto* const token = *tokens;
    if (token->kind == TokenKind::END) {
        throw std::runtime_error(
            "Expected primary expression but the actual is the end of the string");
    }
    const auto c = token->str[0];
    if (c == '(') {
        GoToNextToken(tokens);
        const auto* node = ConstructExpr(tokens, container);
        if ((*tokens)->str != ")") { throw std::runtime_error("Unclosed parenthesis"); }
        GoToNextToken(tokens);
        return node;
    }
    if (token->kind != TokenKind::NUM) {
        throw std::runtime_error(std::string("Expected number expression but the actual is: ") +
                                 token->str.data());
    }
    container.emplace_back(NodeKind::NUM, token, nullptr, nullptr);
    GoToNextToken(tokens);
    return &container.back();
}
const Node* ConstructAST(Token* tokens, std::vector<Node>& container) {
    Token** tokens_ptr = &tokens;
    const auto* root = ConstructExpr(tokens_ptr, container);
    if ((*tokens_ptr)->kind != TokenKind::END) {
        throw std::runtime_error("Finished at not the end of the string (input is not expression)");
    }
    return root;
}
Float Value(const Node* node) {
    using namespace constant::f;
    switch (node->kind) {
        case NodeKind::ADD: return Value(node->l) + Value(node->r);
        case NodeKind::SUB:
            if (node->r == nullptr) {
                return -Value(node->l);
            } else {
                return Value(node->l) - Value(node->r);
            }
        case NodeKind::MUL: return Value(node->l) * Value(node->r);
        case NodeKind::DIV: return Value(node->l) / Value(node->r);
        case NodeKind::NUM:
            if (node->token->str == "pi") {
                return Pi;
            } else {
                return Float(node->token->str);
            }
    }
}
}  // namespace
AST AST::Parse(const std::string& s) {
    auto ast = AST();
    ast.str_ = std::unique_ptr<std::string>(new std::string(s));

    // Tokenize
    ast.tokens_ = std::unique_ptr<std::vector<Token>>(new std::vector<Token>());
    Tokenize(*ast.str_, *ast.tokens_);

    // Construct AST
    ast.nodes_ = std::unique_ptr<std::vector<Node>>(new std::vector<Node>());
    // This `reserve` is important to ensure that the address of node does not change
    ast.nodes_->reserve(ast.tokens_->size());
    ast.root_ = ConstructAST(ast.tokens_->data(), *ast.nodes_);

    return ast;
}
Float AST::Value() const { return ::qrot::Value(root_); }
}  // namespace qrot
