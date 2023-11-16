licenses(["notice"])

cc_library(
    name = "DendAsync",
    hdrs = ["DendAsync.h"],
    deps = [
        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "@gbbs//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "@gbbs//gbbs",
        "//Dendrogram/common:leftist_heap",
        "//Dendrogram/common:skew_heap",
        "//Dendrogram/common:pairing_heap",
        # "//Dendrogram/common:union_find",
    ],
)

cc_library(
    name = "DendSync",
    hdrs = ["DendSync.h"],
    deps = [
        "@gbbs//gbbs",
        "//Dendrogram/common:leftist_heap",
        "//Dendrogram/common:skew_heap",
        "//Dendrogram/common:pairing_heap",
        "//Dendrogram/common:union_find",
    ],
)

cc_binary(
    name = "Dendrogram_main",
    srcs = ["Dendrogram.cc"],
    deps = [
        ":DendAsync",
        ":DendSync",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
