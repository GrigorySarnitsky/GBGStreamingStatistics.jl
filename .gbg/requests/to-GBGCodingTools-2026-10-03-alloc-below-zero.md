Status: open — approved by the owner on 2026-10-03; sent to the GBGCodingTools session the same day.

# To GBGCodingTools: `@alloc_below` cannot assert an exact zero

On Julia 1.12.7 a plain `@allocated merge!(a, b)` on two `Histogram1D`s is 0 in global, function, `try`
and `@testset` scope; inside `@alloc_below 1 merge!(a, b)` the same call reads 48 bytes (96 for the
`Histogram2D` `merge!`, 208 for the keyword form `merge!(c, f; weight = 0.5)`), while `observe!`, `scale!`
and `observe_pairs!` read 0 under the macro. The difference is the macro's own
`@test (@allocated $(esc(expr))) < limit`, whose argument capture allocates.

`RULE.TEST.ALLOC_ZERO` asks a package to pin the true zero where a bound covers another nondeterminism;
through the macro that is not possible. SS therefore pins it with a plain `@test (@allocated …) == 0`
beside each `merge!` row and gives the row a 256-byte bound (`test/allocations.jl`).

Asked: measure inside the macro without the `@test` capture — bind the count to a local first, then test
the local — so that `@alloc_below 1` holds for an operation that allocates nothing.
