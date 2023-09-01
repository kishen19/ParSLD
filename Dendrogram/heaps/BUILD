licenses(["notice"])

cc_library(
    name = "binary_heap",
    hdrs = ["binary_heap.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "leftist_heap",
    hdrs = ["leftist_heap.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "skew_heap",
    hdrs = ["skew_heap.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "pairing_heap",
    hdrs = ["pairing_heap.h"],
    deps = [
        "@gbbs//gbbs",
        ":binary_heap",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
