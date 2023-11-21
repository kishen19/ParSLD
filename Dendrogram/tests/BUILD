licenses(["notice"])

load("@gbbs//internal_tools:build_defs.bzl", "gbbs_cc_test")

gbbs_cc_test(
    name = "test_dendrogram",
    srcs = ["test_dendrogram.cc"],
    deps = [
        "//Dendrogram/SeqUnionFind:Dendrogram",
        "//Dendrogram/ParUnionFind:Dendrogram",
        "//Dendrogram/RCtreeTracing:Dendrogram",
        "//utils:generate_MSF",
        "//utils:caterpillar",
        "//utils:path",
        "//utils:star",
        "//utils:full_binary",
        "//utils:uniform_hook",
        "@gbbs//gbbs:graph",
        "@gbbs//gbbs:gbbs",
        "@gbbs//gbbs:macros",
        "@gbbs//gbbs/unit_tests:graph_test_utils",
        "@googletest//:gtest_main",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
