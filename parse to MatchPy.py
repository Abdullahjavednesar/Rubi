import re

full_string = '''
List[RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Times[9, Power[4, -1]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Times[1, Power[4, -1]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], PosQ[Times[b, d, Power[Times[a, c], -1]]]]]], Plus[Times[-4, Power[Times[5, b, Power[Plus[a, Times[b, x]], Times[5, Power[4, -1]]], Power[Plus[c, Times[d, x]], Times[1, Power[4, -1]]]], -1]], Times[-1, d, Power[Times[5, b], -1], Int[Times[1, Power[Times[Power[Plus[a, Times[b, x]], Times[5, Power[4, -1]]], Power[Plus[c, Times[d, x]], Times[5, Power[4, -1]]]], -1]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Times[5, Power[4, -1]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Times[1, Power[4, -1]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], PosQ[Times[b, d, Power[Times[a, c], -1]]]]]], Plus[Times[-2, Power[Times[b, Power[Plus[a, Times[b, x]], Times[1, Power[4, -1]]], Power[Plus[c, Times[d, x]], Times[1, Power[4, -1]]]], -1]], Times[Plus[Times[b, c], Times[-1, a, d]], Power[Times[2, b], -1], Int[Times[1, Power[Times[Power[Plus[a, Times[b, x]], Times[5, Power[4, -1]]], Power[Plus[c, Times[d, x]], Times[5, Power[4, -1]]]], -1]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Times[3, Power[2, -1]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Times[3, Power[2, -1]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]]]]], Times[x, Power[Times[a, c, Sqrt[Plus[a, Times[b, x]]], Sqrt[Plus[c, Times[d, x]]]], -1]]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Sqrt[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]]], Sqrt[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]]]]], Times[2, Subst[Int[Times[1, Power[Plus[b, Times[-1, d, Power[x, 2]]], -1]], x], x, Times[Sqrt[Plus[a, Times[b, x]]], Power[Sqrt[Plus[c, Times[d, x]]], -1]]]]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Sqrt[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]]], Sqrt[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], PositiveQ[a], ZeroQ[Plus[a, c]]]]], Times[ArcCosh[Times[b, x, Power[a, -1]]], Power[b, -1]]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]]]]], Int[Times[1, Power[Plus[Times[a, c], Times[b, d, Power[x, 2]]], -1]], x]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Times[Plus[Optional[Pattern[a, Blank[]]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Plus[Optional[Pattern[c, Blank[]]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]]], -1]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], NonzeroQ[Plus[Times[b, c], Times[-1, a, d]]]]]], Plus[Times[b, Power[Plus[Times[b, c], Times[-1, a, d]], -1], Int[Times[1, Power[Plus[a, Times[b, x]], -1]], x]], Times[-1, d, Power[Plus[Times[b, c], Times[-1, a, d]], -1], Int[Times[1, Power[Plus[c, Times[d, x]], -1]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], -1]], Pattern[x, Blank[Symbol]]], FreeQ[List[a, b], x]]], Times[Log[RemoveContent[Plus[a, Times[b, x]], x]], Power[b, -1]]], RuleDelayed[HoldPattern[Int[Times[1, Power[Pattern[x, Blank[]], -1]], Pattern[x, Blank[Symbol]]]], Log[x]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], PositiveIntegerQ[Plus[m, Times[1, Power[2, -1]]]]]]], Plus[Times[x, Power[Plus[a, Times[b, x]], m], Power[Plus[c, Times[d, x]], m], Power[Plus[Times[2, m], 1], -1]], Times[2, a, c, m, Power[Plus[Times[2, m], 1], -1], Int[Times[Power[Plus[a, Times[b, x]], Plus[m, -1]], Power[Plus[c, Times[d, x]], Plus[m, -1]]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Pattern[n, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], IntegerQ[Plus[m, Times[1, Power[2, -1]]]], IntegerQ[Plus[n, Times[1, Power[2, -1]]]], Less[m, n, 0]]]], Plus[Times[-1, Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Plus[c, Times[d, x]], Plus[n, 1]], Power[Times[2, a, d, Plus[m, 1]], -1]], Times[Plus[m, n, 2], Power[Times[2, a, Plus[m, 1]], -1], Int[Times[Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Plus[c, Times[d, x]], n]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Pattern[n, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], IntegerQ[Plus[m, Times[1, Power[2, -1]]]], IntegerQ[Plus[n, Times[1, Power[2, -1]]]], Less[0, m, n]]]], Plus[Times[Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Plus[c, Times[d, x]], n], Power[Times[b, Plus[m, n, 1]], -1]], Times[2, c, n, Power[Plus[m, n, 1], -1], Int[Times[Power[Plus[a, Times[b, x]], m], Power[Plus[c, Times[d, x]], Plus[n, -1]]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d, m], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], Not[IntegerQ[Times[2, m]]]]]], Times[Power[Plus[a, Times[b, x]], FracPart[m]], Power[Plus[c, Times[d, x]], FracPart[m]], Power[Power[Plus[Times[a, c], Times[b, d, Power[x, 2]]], FracPart[m]], -1], Int[Power[Plus[Times[a, c], Times[b, d, Power[x, 2]]], m], x]]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], NegativeIntegerQ[Plus[m, Times[3, Power[2, -1]]]]]]], Plus[Times[-1, x, Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Plus[c, Times[d, x]], Plus[m, 1]], Power[Times[2, a, c, Plus[m, 1]], -1]], Times[Plus[Times[2, m], 3], Power[Times[2, a, c, Plus[m, 1]], -1], Int[Times[Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Plus[c, Times[d, x]], Plus[m, 1]]], x]]]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Optional[Pattern[m, Blank[]]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Optional[Pattern[m, Blank[]]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d, m], x], ZeroQ[Plus[Times[b, c], Times[a, d]]], Or[IntegerQ[m], And[PositiveQ[a], PositiveQ[c]]]]]], Int[Power[Plus[Times[a, c], Times[b, d, Power[x, 2]]], m], x]], RuleDelayed[HoldPattern[Condition[Int[Times[Power[Plus[Optional[Pattern[a, Blank[]]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Optional[Pattern[m, Blank[]]]], Power[Plus[Pattern[c, Blank[]], Times[Optional[Pattern[d, Blank[]]], Pattern[x, Blank[]]]], Pattern[n, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, c, d, m, n], x], NonzeroQ[Plus[Times[b, c], Times[-1, a, d]]], ZeroQ[Plus[m, n, 2]], NonzeroQ[Plus[m, 1]]]]], Times[Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Plus[c, Times[d, x]], Plus[n, 1]], Power[Times[Plus[Times[b, c], Times[-1, a, d]], Plus[m, 1]], -1]]], RuleDelayed[HoldPattern[Condition[Int[Power[Plus[Optional[Pattern[a, Blank[]]], Times[Optional[Pattern[b, Blank[]]], Pattern[u, Blank[]]]], Pattern[m, Blank[]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, m], x], LinearQ[u, x], NonzeroQ[Plus[u, Times[-1, x]]]]]], Times[1, Power[Coefficient[u, x, 1], -1], Subst[Int[Power[Plus[a, Times[b, x]], m], x], x, u]]], RuleDelayed[HoldPattern[Condition[Int[Power[Plus[Optional[Pattern[a, Blank[]]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, m], x], NonzeroQ[Plus[m, 1]]]]], Times[Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Times[b, Plus[m, 1]], -1]]], RuleDelayed[HoldPattern[Condition[Int[Power[Pattern[x, Blank[]], Optional[Pattern[m, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[m, x], NonzeroQ[Plus[m, 1]]]]], Times[Power[x, Plus[m, 1]], Power[Plus[m, 1], -1]]]]
'''

full_string = '''
List[RuleDelayed[HoldPattern[Condition[Int[Times[1, Power[Plus[Pattern[a, Blank[]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], -1]], Pattern[x, Blank[Symbol]]], FreeQ[List[a, b], x]]], Times[Log[RemoveContent[Plus[a, Times[b, x]], x]], Power[b, -1]]], RuleDelayed[HoldPattern[Condition[Int[Power[Plus[Optional[Pattern[a, Blank[]]], Times[Optional[Pattern[b, Blank[]]], Pattern[x, Blank[]]]], Pattern[m, Blank[]]], Pattern[x, Blank[Symbol]]], And[FreeQ[List[a, b, m], x], NonzeroQ[Plus[m, 1]]]]], Times[Power[Plus[a, Times[b, x]], Plus[m, 1]], Power[Times[b, Plus[m, 1]], -1]]], RuleDelayed[HoldPattern[Condition[Int[Power[Pattern[x, Blank[]], Optional[Pattern[m, Blank[]]]], Pattern[x, Blank[Symbol]]], And[FreeQ[m, x], NonzeroQ[Plus[m, 1]]]]], Times[Power[x, Plus[m, 1]], Power[Plus[m, 1], -1]]]]
'''

replacements = dict(
        Times="Mul",
        Plus="Add",
        Power="Pow",
)

def parse_full_form(wmexpr):
    out = []
    stack = [out]
    generator = re.finditer(r'[\[\],]', wmexpr)
    last_pos = 0
    for match in generator:
        if match is None:
            break
        position = match.start()
        last_expr = wmexpr[last_pos:position].replace(',', '').replace(']', '').replace('[', '').strip()

        if match.group() == ',':
            if last_expr != '':
                stack[-1].append(last_expr)
        elif match.group() == ']':
            if last_expr != '':
                stack[-1].append(last_expr)
            stack.pop()
            current_pos = stack[-1]
        elif match.group() == '[':
            stack[-1].append([last_expr])
            stack.append(stack[-1][-1])
        last_pos = match.end()
    return out[0]


def generate_sympy_from_parsed(parsed, depth=0):
    out = ""
    if not isinstance(parsed, list):
        return parsed
    if parsed[0] == "If":
        out += "\nif ("+generate_sympy_from_parsed(parsed[1])+"):\n"
        out += shift4(generate_sympy_from_parsed(parsed[2]))
        if len(parsed) > 2:
            out += "\nelse:\n"
            out += shift4(generate_sympy_from_parsed(parsed[3]))
    else:
        if parsed[0] in replacements:
            out += replacements[parsed[0]]
        else:
            out += parsed[0]
        if len(parsed) == 1:
            return out
        out += "("
        out += ", ".join([generate_sympy_from_parsed(i) for i in parsed[1:]])
        out += ")"
    return out


def add_wildcards(string):
    symbols = []

    p = r'(Optional\(Pattern\((\w), Blank\)\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], i[1] + '_')
        symbols.append(i[1])

    p = r'(Pattern\((\w), Blank\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], i[1] + '_')
        symbols.append(i[1])

    p = r'(Pattern\((\w), Blank\(Symbol\)\))'
    matches = re.findall(p, string)
    for i in matches:
        string = string.replace(i[0], i[1] + '_')

    return string, symbols


res = parse_full_form(full_string)
rules = []

for i in res:
    if i[0] == 'RuleDelayed':
        rules.append(i)

parsed = ''

for i in range(0, len(rules)):
    r = rules[i]

    pattern = generate_sympy_from_parsed(r[1][1][1])
    pattern, free_symbols = add_wildcards(pattern)

    condition = generate_sympy_from_parsed(r[1][1][2])
    transformed = generate_sympy_from_parsed(r[2])

    parsed = parsed + 'pattern' + str(i) +' = Pattern(' + pattern + ', cond(' + condition + ', (' + ', '.join(free_symbols) + '))'
    parsed = parsed + '\n' + 'rule' + str(i) +' = ReplacementRule(' + 'pattern1, lambda ' + ', '.join(free_symbols) + ' : ' + transformed + ')\n\n'

print(parsed)
