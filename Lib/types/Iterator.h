
    // struct Iterator{
    //     using iterator_category = std::forward_iterator_tag;
    //     using difference_type   = std::ptrdiff_t;
    //     using value_type        = T;
    //     using pointer           = value_type*;
    //     using reference         = value_type&;
    //     pointer m_ptr;
    //     Iterator(pointer ptr) : m_ptr(ptr) {}
    //     reference operator*() const { return *m_ptr; }
    //     const reference operator*() const { return *m_ptr; }
    //     pointer operator->() { return m_ptr; }
    //     Iterator& operator++() { m_ptr++; return *this; }  
    //     Iterator operator++(int) const{ Iterator tmp = *this; ++(*this); return tmp; }
    //     Iterator& operator--() { m_ptr--; return *this; }  
    //     Iterator operator--(int) { Iterator tmp = *this; --(*this); return tmp; }
    //     constexpr auto operator== (const Iterator& b) const { return m_ptr == b.m_ptr; }
    //     constexpr auto operator!= (const Iterator& b) const { return m_ptr != b.m_ptr; }
    //     constexpr const auto begin() { return v.begin(); }
    //     constexpr const auto end()   { return v.begin()+sz; } 
    //     constexpr const auto begin() const{ return v.begin(); }
    //     constexpr const auto end() const  { return v.begin()+sz; } 
    // };

//https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
