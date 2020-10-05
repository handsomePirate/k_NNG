#pragma once
#include <queue>

template <class T, class S, class C>
const S& Container(const std::priority_queue<T, S, C>& q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
        static const S& Container(const std::priority_queue<T, S, C>& q) {
            return q.*&HackedQueue::c;
        }
    };
    return HackedQueue::Container(q);
}
