#include "Template.hh"

#include "test_main.hh"

TEST(Template, tuple_combinations) {
    using Tuple1 = std::tuple<int, double>;
    using Tuple2 = std::tuple<char, float>;

    using Result = tuple_combinations<Tuple1, Tuple2>;

    using Expected =
        std::tuple<std::pair<int, char>, std::pair<int, float>, std::pair<double, char>, std::pair<double, float> >;
    static_assert(std::is_same<Result, Expected>::value, "Combination result does not match the expected type");
}
