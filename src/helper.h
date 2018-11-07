template <typename T>
void deallocate_container(T& c)
{
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
}
