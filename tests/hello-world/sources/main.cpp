#include "rtk/rtk.hpp"

auto main() -> int {
    scion::log::info("hello", "rtk {}.{}.{}", rtk::version::MAJOR, rtk::version::MINOR, rtk::version::PATCH);

    return 0;
}
